

library(plyr) ## for function ddply
library(gdata) ## for function resample
library(gstat) ## for variogram functions
library(proxy) ## for function dist
library(RandomFields) ## to generate stress landscape
library(tidyverse)
library(magrittr)
library(pracma)
library(dgof)
library(zipfR)
library(hypergeo)
library(orthopolynom)
library(splus2R) ## for function peaks
library(cowplot)

# set printing theme
theme_set(theme_bw())

# Define print colors
red='#DD5144'
green='#1DA462'
blue='#4C8BF5'
yellow='#FFCD46'


do.data = 1
do.analysis = !do.data

## Indicator variable; 1 if running on clusters, 0 if running on personal computer.
hpcc = 0

## Directory to save data to
dirname = '~/DispersalLimitation/Data/20201203_R2'
setwd(dirname)

## Scenarios
## Note: dispersal_scale =  distance within which 95% of a tree's seeds are dispersed, 
## relative to envt "patch size" (range of variogram)
## Note: dispersal_neighborhood: number of sites to the side of an individual of the species with
## avg fecundity encompassing a circle within which all its seeds are dispersed

scenarios =
  crossing(
    age_structure = c(FALSE, TRUE),
    tolerance_steepness = c(0, 15, 100),
    noise_level = c(0, .5),
    mean_fecundity = c(10, 100, 1000),
    dispersal_neighborhood = c(1, 5, 10, Inf),
    replicate = 1:10
  ) %>%
  filter(noise_level == 0 | tolerance_steepness > 0) %>%
  filter(age_structure == FALSE | is.finite(dispersal_neighborhood)) %>%
  mutate(neutrality = tolerance_steepness == 0) %>%
  rowid_to_column(var = 'scenario')

if(do.data){
  
  ## Framing parameters
  {
    ind = 31
    
    RandomField = 
      function(
         L = 50,
         rangepar,
         sillpar,
         nuggetpar,
         seed,
         plot = FALSE
      ) {
        stopifnot(nuggetpar >= 0 & nuggetpar <= 1)
        
        RFoptions(seed = seed)
        
        stress = 
          RFsimulate(
            RMgauss(
              scale = rangepar + 1e-16,
              var = 2 * sillpar * (1 - nuggetpar)
            ) + RMtrend(mean = 0) + RMnugget(var = 2 * sillpar * nuggetpar),
            x = 1:L,
            y = 1:L
          )@data$variable1
        
        if (plot)
          plot(
            raster::raster(matrix(stress, L, L)),
            las = 1,
            xlab = 'x-coordinate',
            ylab = 'y-coordinate'
          )
        
        return(stress)
      
      }
      
    DispersalKernel = function(d, w) exp(-(d / w) ^ 2) / (pi * w ^ 2)
    
    norm_const = function(w){
      n <- 1:1e6
      1 / (1 + 8 * sum(n * exp(-(n / w) ^ 2)))
    } 
    ChessboardDispersal = function(d, w) norm_const(w) * exp(-(d / w) ^ 2)
    CumulativeDispersal = 
      function(d, w){
          n = 0:d
          ChessboardDispersal(0, w) + 8 * sum(n * ChessboardDispersal(n, w))
      } 
    
    Gaussian_variogram = function(dist, range, sill, nugget) 
      (sill - nugget) * (1 - exp(-dist ^ 2 / range ^ 2)) + nugget * (dist > 0)
    
    Gaussian_vgm_optim = function(parms, sample.vgm){
      range = parms[1]
      sill = parms[2]
      nugget = parms[3]
      dist = c(0, sample.vgm$dist)
      observed = c(0, sample.vgm$gamma)
      predicted = (sill - nugget) * (1 - exp(-dist ^ 2 / range ^ 2)) + nugget * (dist > 0)
      sum_squared_errors = sum((observed - predicted) ^ 2)
      return(sum_squared_errors)
    } 
    
    scen = scenarios[ind, ]
    
    list2env(scen, envir = environment())
    
    filename = paste0('scenario_',scenario,'.RData')
    
    set.seed(replicate)
    
    fixed_parameters <-
      tibble(
        landscape_length = 50, ## BCI has 420k trees, i.e. landscape_length = 650
        num_species = 50,   ## BCI has 325 species, https://forestgeo.si.edu/sites/neotropics/barro-colorado-island
        stress_levels = 1 * num_species,
        landscape_autocor_parm = 5,
        delta_t = 1,
        maxtime = 2e5,
        annual_turnover_rate = .1,
        immigration_rate = .01,
        census_interval = 1e3
      )
    
    list2env(fixed_parameters, envir = environment())
    
  }
  
  ## ENVT PARAMETERS: Generate the stress landscape
  {
    field <-
      RandomField(
        L = landscape_length,
        rangepar = landscape_autocor_parm,
        sillpar = 1,
        nuggetpar = noise_level,
        seed = replicate
      )
    
    ## simplify the field by boxing sites into one of stress_levels soil types
    stress_landscape <- 
      1 + findInterval(field, quantile(field, seq(stress_levels) / stress_levels), rightmost.closed = TRUE)
    
    site_info <-
      expand_grid(
        x = seq(landscape_length), 
        y = seq(landscape_length)
      ) %>%
      mutate(stress = stress_landscape) %>% 
      rowid_to_column(var = 'site')
    
  }
  
  ## ENVT PARAMETERS: Calculate envt scale = range of the landscape variogram
  {
    sample.vgm = variogram(stress ~ 1, data = site_info, locations = ~ x + y, width = 1)
    
    fitted.vgm = 
      optim(
        par = c(range = 8, sill = var(stress_landscape), nugget = 1000), 
        fn = Gaussian_vgm_optim,
        sample.vgm = sample.vgm,
        method = "L-BFGS-B",
        lower = c(1, 1, 0)
      )$par
    sill = as.numeric(fitted.vgm[2])		## semivariance at the landscape level --> regional heterogeneity
    nugget = as.numeric(fitted.vgm[3])	## semivariance at distance = 1	--> local uniformity
    range = as.numeric(fitted.vgm[1])		## distance at which the semivariance reaches 63% of the sill
    range95 = sqrt(3) * range   ## distance at which the semivariance reaches 95% of the sill
    
   if(!hpcc){
      plot_variogram <-
        tibble(
          dist = c(0, sample.vgm$dist),
          gamma = c(0, sample.vgm$gamma),
          fitted = Gaussian_variogram(c(0, sample.vgm$dist), range, sill, nugget)
        ) %>% 
        ggplot() + 
        geom_line(aes(dist, fitted), color = blue, size = 1) + 
        geom_point(aes(dist, gamma)) +
        geom_vline(aes(xintercept = range95), color = red) +
        geom_hline(aes(yintercept = sill), color = green) +
        geom_hline(aes(yintercept = nugget), color = yellow) +
        labs(x = 'Distance', y = 'Semivariance') +
        ggtitle('Sample and fitted variogram')
      
      gridExtra::grid.arrange(plot_variogram)
    }
  }
  
  ## SPECIES PARAMETERS: age-specific fertility. 
  ## Trees aged 0 to 5 years have 0 fertility, then between 5 and 10 have 1/2 fertility, then 10 and above have full fertility.
  age_breakdown = c(0, 5, 10)
  if(age_structure == TRUE) fertility_factors = c(0, .5, 1)
  if(age_structure == FALSE) fertility_factors = c(1, 1, 1)
  
  ## SPECIES PARAMETERS: Generate fecundity levels
  {
    fecundity = sort(rlnorm(num_species))
    fecundity = fecundity / mean(fecundity) * mean_fecundity
  }
  
  ## SPECIES PARAMETERS: Calculate tolerance levels
  {
    ## tolerance parameters: 
    ## ensure that the stress threshold (stress at which Tij = 0.5) is the  
    ## highest (lowest) for the most (least) competitive species
    smin = min(stress_landscape)
    smax = max(stress_landscape)
    fmin = min(fecundity)
    fmax = max(fecundity)
    M = matrix(c(smin, smax, fmax, fmin), 2, 2)
    tolpars = solve(M, c(1, 1))
    cs = tolpars[1]
    cf = tolpars[2]
    
    tolerance <- 
      crossing(
        tibble(stress = seq(stress_levels)), 
        tibble(species = seq(num_species), fecundity = fecundity)
      ) %>%
      mutate(tolerance = 1e-16 + .5 * VGAM::erfc(tolerance_steepness * (cs * stress + cf * fecundity - 1)))
    
    s50 <-
      tolerance %>% 
      group_by(species) %>% 
      summarize(s50 = stress[which.min(tolerance - .5)[1]], .groups = 'drop')
    
  }
  
  ## SPECIES PARAMETERS: Reset species parameters for neutral simulation
  {
    if(neutrality){
      fecundity = rep(mean(fecundity), num_species)
      tolerance$fecundity = rep(mean(fecundity), nrow(tolerance))
      tolerance$tolerance = rep(.5, nrow(tolerance))
    }
  }
  
  ## SPECIES PARAMETERS: Define sigma_dispersal based on the envt scale and the parameter "dispersal_scale"
  {
    ## sigma_dispersal: Distance scale used in the dispersal kernel, ensures that all dispersal of the 
    ## avg fec species is within the dispersal_neighborhood
    
    # sigma_dispersal = dispersal_neighborhood / sqrt(log(mean_fecundity))
    
    ## distance_threshold: Distance at which total # of seeds reaching local site for the 
    ## highest-fecundity species < 1: corresponds to a circle within whose radius will fall more than 
    ## (fmax - 1) seeds, such that outside this circle will fall less than 1 seed from the most fecund species.
    ## Calculated via the expression: fmax - fmax * K(d_thresh) < 1, 
    ## where K(d) is the cumulative dispersal, given by K(d) = integrate(2 * pi * r * Kernel(d) * dr, 0, d).
    ## In the case of Gaussian dispersal, we have K(d) = 1 - exp(-(d / sigma_dispersal) ^ 2).
    
    # distance_threshold = sigma_dispersal * sqrt(log(max(fecundity)))
    
    
    ## The parameters below are for the Chessboard (Chebyshev) distance, ie the number of moves
    ## it takes the chess king to get to the square from its original position. It corresponds to 
    ## max(|x - x0|, |y - y0|). For a distance n, there proportion of seeds dispersed within that radius is 
    ## K(0, w) + 8 * sum_1^n(n * K(d, w))
    
    if(dispersal_neighborhood < 1e2){
      sigma_dispersal = 
        uniroot(
          f = function(w){
            n = 1:dispersal_neighborhood
            1 - norm_const(w) * (1 + 8 * sum(n * exp(-(n / w) ^ 2))) - 1 / mean_fecundity
          }, 
          interval = c(1e-3, landscape_length)
        )$root
      
      distance_threshold = 
        uniroot(
          f = function(distance){
            n = 1:distance
            return(1 - norm_const(sigma_dispersal) * (1 + 8 * sum(n * exp(-(n / sigma_dispersal) ^ 2))) - 1 / fmax)
          }, 
          interval = c(1, landscape_length)
        )$root
      
      ## dispersal_scale: the ratio between the distance within which all seeds of the max fec species are dispersed
      ## and the distance where the landscape variogram reaches 95% of its still.
      dispersal_scale = distance_threshold / range95
    }
    
    if(dispersal_neighborhood >= 1e2){
      sigma_dispersal = 0
      
      distance_threshold = Inf
      
      dispersal_scale = Inf
    }
    
  }
  
  ## AUXILIARY PARAMETERS
  {
    dist_matrix = proxy::dist(site_info %>% select(x, y), site_info %>% select(x, y), method = 'maximum')
    
    dispersal_values = sapply(0:(landscape_length - 1), function(d) ChessboardDispersal(d, sigma_dispersal))
    
    disp_matrix = 
      matrix(
        sapply(dist_matrix, function(d) dispersal_values[1+ d]), 
        landscape_length ^ 2, 
        landscape_length ^ 2
      )
    
    tolerance_matrix = 
      unname(
        tolerance %>% 
        select(-fecundity) %>% 
        pivot_wider(names_from = stress, values_from = tolerance) %>% 
        select(-species) %>%
        as.matrix
      )
    
  }
  
  ## Initial conditions
  set.seed(replicate)
  
  {
    species_map = 
      sapply(
        stress_landscape, 
        function(local_stress) 
          gdata::resample(seq(num_species), size = 1, prob = dnorm(local_stress - s50$s50, sd = stress_levels / 10))
      )
    
    data <- 
      bind_cols(
        site_info, 
        species = species_map
      )
    
    current_richness = length(unique(data %>% filter(species > 0) %>% pull(species)))
    
    if(hpcc){
      dat <-
        list(
          scenario = scen,
          stress_landscape = stress_landscape,
          site_info = site_info,
          fixed_parameters = fixed_parameters,
          sample_variogram = sample.vgm,
          fitted_variogram_parameters = 
            tibble(
              range = range, 
              range95 = range95, 
              sill = sill, 
              nugget = nugget
            ),
          spatial_scale_parameters = 
            tibble(
              dispersal_scale = dispersal_scale,
              distance_threshold = distance_threshold,
              sigma_dispersal = sigma_dispersal,
            ),
          distance_threshold = distance_threshold,
          fecundity = tibble(species = seq(num_species), fecundity = fecundity),
          fertility_table = 
            tibble(
              `age threshold` = age_breakdown, 
              `fertility factor` = fertility_factors
            ),
          tolerance = tolerance,
          s50 = s50,
          species_map = tibble(year_0 = species_map)
        )
      save(dat, file = filename)
      rm(dat)
    }
    
    year = 0
    
    ages = rep(10, landscape_length ^ 2)
    
    fec_vec = c(0, fecundity)[1 + data$species] 
    
    fertility = fec_vec * fertility_factors[findInterval(ages, age_breakdown)]
    
    abundances = sapply(seq(num_species), function(species) sum(data$species == species))
  }
  
  ## Dynamics
  set.seed(replicate)
  
  {
    while(year <= maxtime){
      year = year + delta_t
      
      ages = ages + delta_t
      
      if(!hpcc & mod(year, 10) == 0) print(year)
      
      ## death and recruitment
      if(any(data$species > 0)){
        expected_numdeaths = annual_turnover_rate * delta_t * sum(data$species > 0)
        numdeaths = min(2 * expected_numdeaths, rpois(1, lambda = expected_numdeaths))
        
        vacancies = 
          unique(
            c(
              gdata::resample(data$site, size = numdeaths),
              data$site[data$species == 0]
            )
          )
        
        ages[vacancies] = 0
        
        if(dispersal_neighborhood == Inf){
            total_seeds = 
              matrix(
                rpois(
                  num_species * length(vacancies), 
                  lambda = abundances * fecundity / landscape_length ^ 2
                ),
                nrow = num_species
              )
          
            tolerance_vacancies = tolerance_matrix[, stress_landscape[vacancies]]
            
            viable_seeds = 
              matrix(
                ifelse(
                  total_seeds == 0, 
                  0, 
                  rbinom(length(total_seeds), size = total_seeds, prob = as.numeric(tolerance_vacancies))
                ),
                nrow = num_species
              )
            
            recruits = 
              apply(viable_seeds, 2, function(seeds){
                if(sum(seeds) == 0) return(0)
                gdata::resample(seq(num_species), size = 1, prob = seeds)
              })
            
            data$species[vacancies] = recruits
            
            abundances = sapply(seq(num_species), function(species) sum(data$species == species))
            
        }
        
        if(is.finite(dispersal_neighborhood)){
          
            extant_species = intersect(seq(num_species), sort(unique(data$species)))
            occupied_sites = data$site[data$species > 0]
            
            disp_vacancies = disp_matrix[occupied_sites, vacancies]
            
            total_seeds = 
              bind_cols(
                species = data$species[occupied_sites], 
                as_tibble(
                  matrix(
                    rpois(
                      length(disp_vacancies), 
                      lambda = disp_vacancies * fertility[occupied_sites]
                    ), 
                    nrow = length(occupied_sites)
                  )
                )
              ) %>% 
              group_by(species) %>% 
              summarize_all(sum) %>% 
              ungroup %>% 
              select(-species) %>% 
              unlist %>% 
              unname
            
            tolerance_vacancies = 
              tolerance_matrix[
                extant_species, 
                stress_landscape[vacancies]
              ]
            
            viable_seeds = 
              matrix(
                ifelse(
                  total_seeds == 0, 
                  0, 
                  rbinom(length(total_seeds), size = total_seeds, prob = as.numeric(tolerance_vacancies))
                ),
                nrow = length(extant_species)
              )
            
            recruits = 
              apply(viable_seeds, 2, function(seeds){
                if(sum(seeds) == 0) return(0)
                sample(extant_species, size = 1, prob = seeds)
              })
            
            data$species[vacancies] = recruits
            
            fec_vec[vacancies] = c(0, fecundity)[1 + recruits]
            
        }
        
      }
      
      ## immigration
      if(runif(1) <= immigration_rate & any(data$species == 0)){
        spot = gdata::resample(data$site[data$species == 0], size = 1)
        immigrant_species = 
          gdata::resample(
            seq(num_species), 
            size = 1,
            prob = fecundity * tolerance_matrix[, data$stress[spot]]
          )
        data$species[spot] = immigrant_species
        fec_vec[spot] = fecundity[immigrant_species]
      } 
      
      extant_species = intersect(seq(num_species), sort(unique(data$species)))
      current_richness = length(extant_species)
      
      fertility = fec_vec * fertility_factors[findInterval(ages, age_breakdown)]
      
      if(
        mod(
          year, 
          ifelse(
            hpcc, 
            census_interval, 
            ifelse(is.finite(dispersal_neighborhood), 10, 1000)
           )
        ) == 0 | 
        current_richness == 0
        ){
        if(!hpcc){
          plot1 =
            data %>% 
            ggplot(aes(x, y, fill = species)) + 
            geom_tile() +
            ggtitle(
              paste(
                'Richness =', 
                current_richness,
                '      Empty sites =',
                sum(data$species == 0)
              )
            ) +
            theme(aspect.ratio = 1)
          plot2 = 
            data %>% 
            ggplot(aes(stress, species)) + 
            geom_point() +
            ggtitle('Species to stress mapping') +
            theme(aspect.ratio = 1)
          
          gridExtra::grid.arrange(plot1, plot2, nrow = 1) 
        }
        if(hpcc){
          dat = get(load(filename))
          dat$ages = ages
          dat$species_map = bind_cols(dat$species_map, dummy_name = data$species)
          names(dat$species_map)[ncol(dat$species_map)] = paste0('year_', year)
          save(dat, file = filename)
          rm(dat)
        }
      } 
      
    }  
  }
  
}


if(do.analysis){
  
  read_data = 0
  Fig_landscape = 0
  Fig_ricnhess = 0
  Fig_tolerance = 0
  Fig_species_stress = 0
  Fig_stationarity = 1
  Figs_age_structure = 0
  Fig_niche_overlap = 0
  
  ## Read simulation data
  if(read_data){
   filenames = paste0('scenario_', 1:1050,'.RData')
   
   results <-
      ddply(tibble(filename = filenames), .(filename), function(df){
        filename = df$filename
        dat = get(load(filename))
        species_map = as.numeric(unlist(dat$species_map[, ncol(dat$species_map)]))
        empty_sites = sum(species_map == 0)
        extant_species = species_map[species_map > 0]
        richness = length(unique(extant_species))
        dat$scenario$richness = richness
        dat$scenario$years = as.numeric(unlist(strsplit(rev(names(dat$species_map))[1], 'year_'))[2])
        dat$scenario$empty_sites = empty_sites
        dat$scenario = cbind(dat$scenario, dat$spatial_scale_parameters)
        return(dat$scenario)
      }) %>%
      as_tibble %>%
      arrange(scenario) %>%
      select(-filename)
    
    if('dispersal_neighborhood' %in% names(results)){
      results = 
        results %>% 
        mutate(dispersal_scale = dispersal_neighborhood)
    } 
    
    
    results_summary <-
      results %>%
      filter(years > 0 & dispersal_scale > .1) %>%
      group_by(tolerance_steepness, noise_level, mean_fecundity, dispersal_scale, age_structure) %>%
      summarize_at(
        c('richness', 'empty_sites', 'years'),
        .funs = list(mean = mean, sd = sd)
      ) %>%
      ungroup
    
  }
  
  ## Figure - Tolerance
  if(Fig_tolerance){
    dat0 = get(load(filenames[1]))    ## neutrality
    dat15 = get(load(filenames[10]))  ## smooth tolerance transition
    dat100 = get(load(filenames[19])) ## sharp tolerance transition
    
    tol = 
      bind_rows(
        dat0$tolerance %>% mutate(tolerance_steepness = 0),
        dat15$tolerance %>% mutate(tolerance_steepness = 15),
        dat100$tolerance %>% mutate(tolerance_steepness = 100)
      )
    
    main_dat15 =
      tol %>%
      filter(tolerance_steepness == 15) %>%
      ggplot(aes(species, stress, fill = tolerance)) + 
      geom_tile() +
      theme(aspect.ratio = 1) +
      ggtitle('A') +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    inset_dat15 =
      tol %>%
      filter(tolerance_steepness == 15 & species == 43) %>%
      ggplot(aes(stress, tolerance)) +
      geom_line() +
      geom_point() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    main_dat100 =
      tol %>%
      filter(tolerance_steepness == 100) %>%
      ggplot(aes(species, stress, fill = tolerance)) + 
      geom_tile() +
      theme(aspect.ratio = 1) +
      ggtitle('B') +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    inset_dat100 =
      tol %>%
      filter(tolerance_steepness == 100 & species == 43) %>%
      ggplot(aes(stress, tolerance)) +
      geom_line() +
      geom_point() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    
    fig1a = 
      ggdraw() +
      draw_plot(main_dat15) +
      draw_plot(inset_dat15, x = .15, y = .2, width = .35, height = .3)
    
    fig1b = 
      ggdraw() +
      draw_plot(main_dat100) +
      draw_plot(inset_dat100, x = .15, y = .2, width = .35, height = .3)
    
    gridExtra::grid.arrange(fig1a, fig1b, nrow = 1)
  }
  
  ## Figure - Landscape
  if(Fig_landscape){
    
    ind =
      scenarios %>% 
      filter(
        replicate == 3 & 
        tolerance_steepness == 100 & 
        mean_fecundity == 100 & 
        dispersal_neighborhood == 5
      ) %>% 
      pull(scenario)
    
    clean = ind[1]
    noisy = ind[2]
    
    Gaussian_variogram = function(dist, range, sill, nugget) 
      (sill - nugget) * (1 - exp(-dist ^ 2 / range ^ 2)) + nugget * (dist > 0)
    
    dat_clean = get(load(filenames[clean]))
    L = dat_clean$fixed_parameters$landscape_length
    
    plot_landscape_clean = 
      ggimage(matrix(dat_clean$stress_landscape, L, L)) +
      theme(aspect.ratio = 1) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      scale_fill_gradient2(midpoint = 25, low = "darkgreen", mid = 'yellow', high = 'brown') +
      labs(
        x = 'x coordinate',
        y = 'y coordinate',
        fill = 'stress'
      ) +
      ggtitle('A')
    
    plot_variogram_clean <-
      bind_rows(
        tibble(np = NA, dist = 0, gamma = 0, dir.hor = 0, dir.ver = 0, id = 'var1'),
        dat_clean$sample_variogram
      ) %>%
      mutate(
        fitted = 
          Gaussian_variogram(
            c(dist), 
            range = dat_clean$fitted_variogram_parameters$range, 
            sill = dat_clean$fitted_variogram_parameters$sill, 
            nugget = dat_clean$fitted_variogram_parameters$nugget
            )
        ) %>%
      ggplot() + 
      geom_line(aes(dist, fitted), color = blue, size = 1) + 
      geom_point(aes(dist, gamma)) +
      geom_vline(aes(xintercept = dat_clean$fitted_variogram_parameters$range95), color = red) +
      geom_hline(aes(yintercept = dat_clean$fitted_variogram_parameters$sill), color = green) +
      geom_hline(aes(yintercept = dat_clean$fitted_variogram_parameters$nugget), color = yellow) +
      labs(x = 'distance', y = 'semivariance') +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(aspect.ratio = 1) +
      ggtitle('B')
    
    dat_noisy = get(load(filenames[noisy]))
    L = dat_noisy$fixed_parameters$landscape_length
    
    plot_landscape_noisy = 
      ggimage(matrix(dat_noisy$stress_landscape, L, L)) +
      theme(aspect.ratio = 1) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      scale_fill_gradient2(midpoint = 25, low = "darkgreen", mid = 'yellow', high = 'brown') +
      labs(
        x = 'x coordinate',
        y = 'y coordinate',
        fill = 'stress'
      ) +
      ggtitle('C')
    
    plot_variogram_noisy <-
      bind_rows(
        tibble(np = NA, dist = 0, gamma = 0, dir.hor = 0, dir.ver = 0, id = 'var1'),
        dat_noisy$sample_variogram
      ) %>%
      mutate(
        fitted = 
          Gaussian_variogram(
            c(dist), 
            range = dat_noisy$fitted_variogram_parameters$range, 
            sill = dat_noisy$fitted_variogram_parameters$sill, 
            nugget = dat_noisy$fitted_variogram_parameters$nugget
          )
      ) %>%
      ggplot() + 
      geom_line(aes(dist, fitted), color = blue, size = 1) + 
      geom_point(aes(dist, gamma)) +
      geom_vline(aes(xintercept = dat_noisy$fitted_variogram_parameters$range95), color = red) +
      geom_hline(aes(yintercept = dat_noisy$fitted_variogram_parameters$sill), color = green) +
      geom_hline(aes(yintercept = dat_noisy$fitted_variogram_parameters$nugget), color = yellow) +
      labs(x = 'distance', y = 'semivariance') +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(aspect.ratio = 1) +
      ggtitle('D')
    
    
    gridExtra::grid.arrange(
      plot_landscape_clean, 
      plot_variogram_clean,
      plot_landscape_noisy,
      plot_variogram_noisy,
      nrow = 2, 
      widths = c(1, .875))
  }
  
  ## Figure - Richness
  if(Fig_ricnhess){
    
    {
      
      average_variogram_range95 = 8.961648
      
      noise_labs = c('smooth patches', 'noisy patches')
      
      data_for_plotting = 
        results_summary %>%
        bind_rows(
          results_summary %>%
            filter(tolerance_steepness == 0) %>%
            mutate(noise_level = .5)
        ) %>%
        filter(age_structure == FALSE) %>%
        mutate(
          dispersal_scale = factor(round(dispersal_scale / average_variogram_range95, 1)),
          tolerance_steepness = factor(tolerance_steepness),
          noise_level = noise_labs[match(noise_level, c(0, .5))]
        ) %>%
        mutate(noise_level = factor(noise_level, levels = noise_labs))
      
      plot_richness <-
        data_for_plotting %>%
        ggplot(
          aes(
            x = dispersal_scale,
            y = richness_mean, 
            fill = tolerance_steepness,
            group = tolerance_steepness
          )
        ) +
        geom_col(position = 'dodge', color = 'black') +
        geom_errorbar(
          aes(
            x = dispersal_scale, 
            ymin = richness_mean - richness_sd, 
            ymax = richness_mean + richness_sd
          ),
          position = position_dodge(width = .85),
          width = .2
        ) +
        facet_grid(
          noise_level ~ mean_fecundity,
          labeller = 
            label_bquote(
              cols = paste(fec[mean], ' = ',  .(mean_fecundity))
            )
        ) +
        theme(strip.background = element_rect(fill = yellow)) +
        labs(x = 'dispersal scale / range of environmental variation', y = 'richness') +
        scale_fill_discrete(name = 'tolerance \n gradient', labels = c('neutral', 'gradual', 'steep')) +
        theme(aspect.ratio = 1)
     
      plot_emptysites <-
        results_summary %>%
        bind_rows(
          results_summary %>%
            filter(tolerance_steepness == 0) %>%
            mutate(noise_level = 0.5)
        ) %>%
        filter(age_structure == FALSE) %>%
        mutate(
          dispersal_scale = factor(round(dispersal_scale / average_variogram_range95, 1)),
          tolerance_steepness = c('neutral', 'gradual', 'steep')[match(tolerance_steepness, c(0, 15, 100))],
          empty_sites_mean = 100 * empty_sites_mean / 2500,
          empty_sites_sd = 100 * empty_sites_sd / 2500,
          noise_level = noise_labs[match(noise_level, c(0, .5))]
        ) %>%
        mutate(
          tolerance_steepness = factor(tolerance_steepness, levels = c('neutral', 'gradual', 'steep')),
          noise_level = factor(noise_level, levels = noise_labs)
        ) %>%
        ggplot(
          aes(
            x = dispersal_scale, 
            y = empty_sites_mean,
            fill = tolerance_steepness,
            group = tolerance_steepness
          )
        ) +
        geom_col(position = 'dodge', color = 'black') +
        geom_errorbar(
          aes(
            x = dispersal_scale, 
            ymin = empty_sites_mean - empty_sites_sd, 
            ymax = empty_sites_mean + empty_sites_sd
          ),
          position = position_dodge(width = .85),
          width = .2
        ) +
        facet_grid(
          noise_level ~ mean_fecundity, 
          labeller = 
            label_bquote(
              cols = paste(fec[mean], ' = ',  .(mean_fecundity))
            )
        ) +
        theme(strip.background = element_rect(fill = yellow)) +
        labs(
          x = 'dispersal scale / range of environmental variation',
          y = 'empty sites (%)',
          fill = 'tolerance \n gradient'
        ) +
        theme(aspect.ratio = 1)
      
      
    }
    
    plot_agestructure_ratio =
      results_summary %>%
      filter(tolerance_steepness > 0) %>%
      select(
        tolerance_steepness, 
        noise_level, 
        mean_fecundity, 
        dispersal_scale, 
        age_structure, 
        richness_mean,
        richness_sd
      ) %>%
      mutate(age_structure = ifelse(age_structure, 'true', 'false')) %>%
      pivot_wider(names_from = age_structure, values_from = c(richness_mean, richness_sd)) %>%
      filter(is.finite(dispersal_scale)) %>%
      mutate(
        rich_ratio = 100 * (richness_mean_true / richness_mean_false - 1),
        rich_ratio_se = rich_ratio * 
          sqrt((richness_sd_false / richness_mean_false) ^ 2 + (richness_sd_true / richness_mean_true) ^ 2)
      ) %>%
      mutate(
        dispersal_scale = factor(round(dispersal_scale / average_variogram_range95, 1)),
        tolerance_steepness = factor(tolerance_steepness),
        mean_fecundity = factor(mean_fecundity),
        noise_level = c('smooth', 'noisy')[match(noise_level, c(0, .5))]
      ) %>%
      mutate(noise_level = factor(noise_level, levels = c('smooth', 'noisy'))) %>%
      ggplot(
        aes(
          x = dispersal_scale,
          y = rich_ratio, 
          fill = tolerance_steepness,
          group = tolerance_steepness
        )
      ) +
      geom_col(position = 'dodge', color = 'black') +
      geom_errorbar(
        aes(
          x = dispersal_scale,
          ymin = rich_ratio - rich_ratio_se,
          ymax = rich_ratio + rich_ratio_se
        ),
        position = position_dodge(width = .85),
        width = .2
      ) +
      facet_grid(noise_level ~ mean_fecundity) +
      theme(strip.background = element_rect(fill = yellow)) +
      labs(x = 'dispersal scale / range of environmental variation', y = 'richness (% difference)') +
      scale_fill_discrete(name = 'tolerance \n gradient', labels = c('gradual', 'steep')) +
      theme(aspect.ratio = 1)
    
    
  }
  
  ## Figure - Species_Stress Association
  if(Fig_species_stress){
    ind = 
      scenarios %>% 
      filter(
        replicate == 1 & 
          dispersal_neighborhood == 1 & 
          mean_fecundity == 100 & 
          noise_level == 0
      ) %>% 
      pull(scenario)
    
    dat0 = get(load(filenames[ind[1]]))
    dat15 = get(load(filenames[ind[2]]))
    dat100 = get(load(filenames[ind[3]]))
    
    the_regime = 'sharp'
    
    chosen_dat = if(the_regime == 'sharp') dat100 else dat15
    
    df =
      dat0$site_info %>%
      add_column(neutral = unlist(dat0$species_map[, ncol(dat0$species_map)])) %>%
      add_column(smooth = unlist(dat15$species_map[, ncol(dat15$species_map)])) %>%
      add_column(sharp = unlist(dat100$species_map[, ncol(dat100$species_map)])) %>%
      pivot_longer(
        cols = c(neutral, smooth, sharp), 
        names_to = 'regime', 
        values_to = 'species'
      ) %>%
      mutate(species = na_if(species, y = 0))
    
    plot_species_map = 
      df %>%
      filter(regime == the_regime) %>% 
      ggplot(aes(x, y, fill = 50 - species)) +
      geom_tile() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      scale_fill_gradient(low = 'black', high = 'deepskyblue') +
      theme(aspect.ratio = 1) +
      labs(
        x = 'x coordinate',
        y = 'y coordinate',
        fill = 'tolerance \n index'
      )
    
    plot_correlation =
      df %>%
      filter(regime == the_regime) %>% 
      group_by(stress, species) %>% 
      count %>% 
      ungroup %>% 
      ggplot(aes(species, stress)) + 
      geom_point() +
      theme(aspect.ratio = 1)
    
    
    plot_association =
      df %>% 
      left_join(
        dat100$tolerance %>% 
          select(species, fecundity), 
        by = 'species'
      ) %>%
      filter(regime == the_regime) %>% 
      group_by(fecundity) %>%
      summarize(mean_stress = mean(stress), .groups = 'drop') %>%
      ggplot(aes(fecundity, mean_stress)) +
      geom_point() +
      theme(aspect.ratio = 1) +
      labs(y = 'mean stress at residency') # +
    # scale_x_log10()
    
    foo = 
      df %>% 
      filter(regime == the_regime) %>% 
      group_by(species) %>%
      summarize(mean_stress = mean(stress), .groups = 'drop') %>%
      left_join(chosen_dat$s50, by = 'species')
    
    mod = lm(mean_stress ~ s50, data = foo %>% filter(s50 < 50))
    print(mod)
    
    plot_meanstress_s50 =
      foo %>%
      filter(s50 < 50) %>%
      ggplot(aes(s50, mean_stress)) +
      geom_smooth(method = 'lm', color = rgb(153 / 255, 0, 0)) +
      geom_point() +
      labs(
        x = 'stress at 50% tolerance',
        y = 'mean stress at residency'
      ) +
      theme(aspect.ratio = 1)
    
    gridExtra::grid.arrange(
      plot_species_map + ggtitle('A'),
      # plot_correlation + ggtitle('B'),
      plot_association + ggtitle('B'),
      plot_meanstress_s50 + ggtitle('C'),
      nrow = 1,
      widths = c(1, .775, .775)
    )
    
    do_SI_fig = 1
    
    if(do_SI_fig){
      
      reread_data = FALSE
      if(reread_data){
        bar = NULL
        for(char in filenames[1:1050]){
          dat = get(load(char))
          species = unlist(dat$species_map[,ncol(dat$species_map)])
          extant_species = species[species > 0]
          bar =
            rbind(
              bar,
              dat$scenario %>%
                bind_cols(
                  dat$site_info %>%
                    add_column(species = species)
                ) %>%
                left_join(dat$s50 %>% select(-.groups), by = 'species')
            )
        } 
      } else{
        bar = get(load('final_state_all_scenarios.rdata'))
      }
      
      
      bar_sum =
        bar %>%
        filter(species > 0) %>%
        group_by(
          scenario, 
          age_structure,
          noise_level,
          tolerance_steepness, 
          dispersal_neighborhood, 
          mean_fecundity, 
          species, 
          s50
        ) %>%
        summarize(stress = mean(stress), .groups = 'drop')
      
      plotting_dat = 
        bar_sum %>% 
        filter(
          tolerance_steepness == 100,
          s50 < 50, s50 > 10,
          mean_fecundity == 100,
          age_structure == FALSE,
          noise_level == 0
        ) 
      
      plot_meanstress_species_all = 
        plotting_dat %>% 
        mutate(dispersal_neighborhood = factor(dispersal_neighborhood)) %>%
        ggplot(
          aes(
            species, 
            stress, 
            group = dispersal_neighborhood, 
            color = dispersal_neighborhood
          )
        ) + 
        geom_point() + 
        geom_smooth(se = FALSE, method = 'loess') +
        ylab('mean stress at residency') +
        theme(aspect.ratio = 1) +
        theme(legend.position = 'none')
      
      average_variogram_range95 = 8.961648
      
      plot_meanstress_s50_all = 
        plotting_dat %>%
        mutate(dispersal_scale = factor(round(dispersal_neighborhood / average_variogram_range95, 1))) %>%
        ggplot(
          aes(
            s50, 
            stress, 
            group = dispersal_scale, 
            color = dispersal_scale
          )
        ) + 
        geom_point() + 
        geom_smooth(se = FALSE, method = 'lm') +
        labs(
          x = 'stress at 50% tolerance',
          y = 'mean stress at residency',
          color = 'disp / envt scale'
        ) +
        theme(aspect.ratio = 1)
      
      gridExtra::grid.arrange(
        plot_meanstress_species_all + ggtitle('A'),
        plot_meanstress_s50_all + ggtitle('B'),
        nrow = 1,
        widths = c(.755, 1)
      )
      
      slope_intercept =
        plotting_dat %>%
        group_by(dispersal_neighborhood) %>%
        summarize(
          intercept = coef(lm(stress ~ s50))[1],
          slope = coef(lm(stress ~ s50))[2],
          sslope = coef(summary(lm(stress ~ s50)))[2, 2],
          .groups = 'drop'
        )
      
      print(slope_intercept)
      
      data_for_plotting =
        bar %>% 
        filter(
          noise_level == 0 & 
            tolerance_steepness == 100 & 
            mean_fecundity == 100 & 
            replicate == 1
        )
      
      disp.labels = c('0.1', '0.6', '1.1', 'Inf') 
      names(disp.labels) = c('1', '5', '10', 'Inf')
      
      plot_association_SI = 
        data_for_plotting %>%
        ggplot(aes(x, y, fill = 50 - species)) +
        geom_tile() +
        facet_wrap(~dispersal_neighborhood, nrow = 2, labeller = labeller(dispersal_neighborhood = disp.labels)) +
        theme(aspect.ratio = 1) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        scale_fill_gradient(low = 'black', high = 'deepskyblue') +
        labs(
          x = 'x coordinate',
          y = 'y coordinate',
          fill = 'tolerance \n index'
        ) +
        theme(strip.background = element_rect(fill = "#FFCD46")) +
        ggtitle('B')
      
      plot_landscape_SI = 
        data_for_plotting %>%
        filter(dispersal_neighborhood == 1) %>%
        ggplot(aes(x, y, fill = stress)) +
        geom_tile() + 
        theme(aspect.ratio = 1) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        scale_fill_gradient2(midpoint = 25, low = "darkgreen", mid = 'yellow', high = 'brown') +
        ggtitle('A') +
        labs(
          x = 'x coordinate',
          y = 'y coordinate'
        )
      
      
      print_richness =
        data_for_plotting %>%
        group_by(dispersal_neighborhood) %>%
        summarize(richness = length(unique(species)), .groups = 'drop')
      
      gridExtra::grid.arrange(plot_landscape_SI, plot_association_SI, nrow = 1)
      
      egg::ggarrange(
        plot_landscape_SI, 
        plot_association_SI, 
        plot_association + ggtitle('C'),
        plot_meanstress_s50 + ggtitle('D'),
        nrow = 2
      )
      
    }
    
    
  }
  
  ## Figure - Stationarity
  if(Fig_stationarity){
    setwd('~/DispersalLimitation/Data/20201203_R2/')
    filenames = paste0('scenario_', 1:1050, '.rdata')
    
    average_variogram_range95 = 8.961648
    
    ts_tibble = NULL
    for(char in filenames){
      dat = get(load(char))
      richness_vector = 
        dat$species_map %>%
        as.matrix %>%
        apply(., 2, function(v) length(unique(v[v > 0]))) %>% 
        as.numeric
      
      timeseries = 
        dat$scenario %>%
        bind_cols(
          year = seq(0, 1000 * (length(richness_vector) - 1), by = 1000),
          richness = richness_vector
        )
      
      ts_tibble = rbind(ts_tibble, timeseries)
      
    }
    
    trend_pvalue =
      ts_tibble %>% 
      filter(year > 100e3) %>%
      group_by(year, mean_fecundity, dispersal_neighborhood, tolerance_steepness) %>%
      summarize(richness = mean(richness), .groups = 'drop') %>%
      group_by(mean_fecundity, dispersal_neighborhood, tolerance_steepness) %>%
      summarize(
        pval = 
          richness %>% 
          snpar::cs.test() %$% 
          p.value,
        .groups = 'drop'
      )
    
    res = 
      ts_tibble %>%
      select(mean_fecundity, dispersal_neighborhood, tolerance_steepness) %>%
      unique %>%
      left_join(trend_pvalue) %>%
      mutate(stationary = 1 * (pval > .05))
    
    res %>% 
      group_by(tolerance_steepness) %>%
      summarize(
        cases = sum(is.finite(pval)),
        significance = mean(pval < .05, na.rm = TRUE),
        .groups = 'drop')
    
    res %>%
      select(-pval) %>%
      pivot_wider(names_from = mean_fecundity, values_from = stationary) %>%
      arrange(tolerance_steepness)
    
    tol_labs = c('tol. grad. = neutral', 'tol. grad. = gradual', 'tol. grad. = steep')
    
    plotting_data =
      ts_tibble %>%
      filter(
        age_structure == FALSE, 
        year >= 0e3, 
        noise_level == 0
      ) %>%
      mutate(
        tolerance_steepness = tol_labs[match(tolerance_steepness, c(0, 15, 100))],
        tolerance_steepness = factor(tolerance_steepness, levels = tol_labs),
        noise_level = as.factor(noise_level),
        dispersal_scale = factor(round(dispersal_neighborhood / average_variogram_range95, 1))
      ) %>%
      group_by(year, mean_fecundity, dispersal_scale, tolerance_steepness) %>%
      summarize(richness = mean(richness), .groups = 'drop') %>%
      group_by(mean_fecundity, dispersal_scale, tolerance_steepness) %>%
      mutate(loss_rate = c(diff(richness), 0)) %>% 
      ungroup
    
    
    plot_timeseries = 
      plotting_data %>%
      ggplot(aes(year / 1e3, richness, color = dispersal_scale)) +
      geom_line(size = 1) +
      facet_grid(
        tolerance_steepness ~ mean_fecundity, 
        labeller = 
          label_bquote(
            cols = paste(fec[mean], ' = ',  .(mean_fecundity))
          )
      ) +
      theme(strip.background = element_rect(fill = "#FFCD46")) +
      labs(
        x = 'year [thousands]', 
        color = 'dispersal \n scale / \n range of \n environmental \n variation'
      ) +
      theme(aspect.ratio = 1) +
      expand_limits(y = 0) 
      # scale_x_log10() +
      # scale_y_log10()
    
    plot_timeseries_lossrate = 
      plotting_data %>%
      filter(year > 1e3) %>%
      ggplot(aes(year / 1e3, loss_rate, color = dispersal_scale)) +
      geom_line(size = 1) +
      facet_grid(
        tolerance_steepness ~ mean_fecundity, 
        labeller = 
          label_bquote(
            cols = paste(fec[mean], ' = ',  .(mean_fecundity))
          )
      ) +
      theme(strip.background = element_rect(fill = "#FFCD46")) +
      labs(
        x = 'year (thousands)', 
        color = 'dispersal \n scale / \n range of \n environmental \n variation'
      ) +
      theme(aspect.ratio = 1) +
      expand_limits(y = 0) +
      scale_x_log10() +
      labs(
        x = 'year [thousands]',
        y = 'loss rate [species per thousand years]'
      )
    
    
    gridExtra::grid.arrange(plot_timeseries)
    
    table = 
      plotting_data %>%
      group_by(tolerance_steepness, mean_fecundity, dispersal_scale) %>%
      mutate(year = year / 100e3) %>%
      summarize(
        slope = round(coef(lm(richness ~ year))[2], 2),
        .groups = 'drop'
      ) %>%
      pivot_wider(names_from = mean_fecundity, values_from = slope)
    
    xtable::xtable(table)
    
  }
  
  ## Figures - Age strucure
  if(Figs_age_structure){
    
    library(caret)
    library(randomForest)
    
    ## Plot empirical age distribution for scenario 23
    if(FALSE){
      
      ind = 
        scenarios %>% 
        filter(
          age_structure == TRUE,
          tolerance_steepness == 100,
          noise_level == 0,
          mean_fecundity == 100,
          dispersal_neighborhood == 5,
          replicate == 1
        ) %$%
        scenario
      
      dat = get(load(paste0('~/DispersalLimitation/Data/20201203_R2/scenario_', ind,'.RData')))
      
      plot_ages = 
        tibble(age = dat$ages) %>% 
        arrange(age) %>% 
        rowid_to_column(var = 'cdf') %>%
        ggplot(aes(age, cdf / max(cdf))) +
        geom_line() +
        geom_vline(aes(xintercept = 0), color = rgb(153 / 255, 0, 0)) +
        geom_vline(aes(xintercept = 5), color = rgb(153 / 255, 0, 0)) +
        geom_vline(aes(xintercept = 10), color = rgb(153 / 255, 0, 0)) +
        ylab('cumulative distribution') +
        theme(aspect.ratio = 1) +
        scale_x_log10()
      
      plot_ages
    }
    
    ## Predict richness using unstructured data
    if(FALSE){
      setwd('~/DispersalLimitation/Data/20201203_R2/')
      
      bundle_unstructured = NULL
      for(filename in paste0('scenario_', 1:600, '.rdata')){
        dat = get(load(filename))
        
        x = 
          tibble(
            site = 1:(dat$fixed_parameters$landscape_length ^ 2),
            species = as.numeric(unlist(dat$species_map[, ncol(dat$species_map)])),
            age = dat$ages,
            fertility = 1,
            fecundity = c(0, dat$fecundity$fecundity)[1 + species]
          )
        
        y = 
          dat$scenario %>%
          bind_cols(
            x %>%
              summarize(
                nominal_mean_f = dat$scenario$mean_fecundity,
                unweighted_mean_f = mean(fecundity),
                weighted_mean_f = mean(fecundity * fertility),
                richness = length(setdiff(unique(species), 0))
              )
          )
        
        bundle_unstructured = rbind(bundle_unstructured, y)
        
      }
      
      bundle_structured = NULL
      for(filename in paste0('scenario_', 601:1050, '.rdata')){
        dat = get(load(filename))
        
        x = 
          tibble(
            site = 1:(dat$fixed_parameters$landscape_length ^ 2),
            species = as.numeric(unlist(dat$species_map[, ncol(dat$species_map)])),
            age = dat$ages,
            fertility = dat$fertility_table$`fertility factor`[findInterval(age, dat$fertility_table$`age threshold`)],
            fecundity = c(0, dat$fecundity$fecundity)[1 + species]
          )
        
        y = 
          dat$scenario %>%
          bind_cols(
            x %>%
              summarize(
                nominal_mean_f = dat$scenario$mean_fecundity,
                unweighted_mean_f = mean(fecundity),
                weighted_mean_f = mean(fecundity * fertility),
                richness = length(setdiff(unique(species), 0))
              )
          )
        
        bundle_structured = rbind(bundle_structured, y)
        
      }
      
      
      data_richness = 
        bundle_structured %>%
        mutate(age_structure = 'structured') %>%
        bind_rows(
          bundle_unstructured %>%
            mutate(age_structure = 'unstructured')
        ) %>%
        filter(neutrality == FALSE & is.finite(dispersal_neighborhood)) %>%
        mutate(
          tolerance_steepness = as.factor(tolerance_steepness),
          noise_level = as.factor(noise_level),
          seed_output = log10(weighted_mean_f)
        )
      
      train_data = 
        data_richness %>% 
        filter(age_structure == 'unstructured')
      
      new_data = 
        data_richness %>% 
        filter(age_structure == 'structured')
      
      
      ## Linear regression of log10(weighted mean fecundity) ~ log10(nominal mean fecundity)
      regression = 
        lm(
          log10(weighted_mean_f) ~ log10(nominal_mean_f), 
          data = new_data
        )
      
      summary(regression)
      
      plot_correspondence = 
        new_data %>%
        ggplot(aes(log10(unweighted_mean_f), log10(weighted_mean_f))) + 
        geom_point() + 
        geom_smooth(method = 'lm') +
        geom_abline(aes(slope = 1, intercept = 0))
      
      gridExtra::grid.arrange(plot_correspondence)
      
      
      ## Variable importance
      mod_rf_unstructured = 
        randomForest(
          richness ~ 
            tolerance_steepness + 
            noise_level + 
            seed_output + 
            dispersal_neighborhood, 
          importance = TRUE,
          data = train_data
        )
      
      mod_glm_unstructured = 
        glm(
          richness ~ 
            tolerance_steepness + 
            noise_level + 
            seed_output + 
            dispersal_neighborhood, 
          family = gaussian,
          data = train_data
        )
      
      varimp_rf_unstructured = varImp(mod_rf_unstructured)
      writeLines('Unstructured Random Forest')
      varimp_rf_unstructured %>% arrange(desc(Overall))
      
      varimp_glm_unstructured = varImp(mod_glm_unstructured)
      writeLines('Unstructured GLM')
      varimp_glm_unstructured %>% arrange(desc(Overall))
      
      
      mod_rf_structured = 
        randomForest(
          richness ~ 
            tolerance_steepness + 
            noise_level + 
            seed_output + 
            dispersal_neighborhood, 
          importance = TRUE,
          data = new_data
        )
      
      mod_glm_structured = 
        glm(
          richness ~ 
            tolerance_steepness + 
            noise_level + 
            seed_output + 
            dispersal_neighborhood, 
          family = gaussian,
          data = new_data
        )
      
      varimp_rf_structured = varImp(mod_rf_structured)
      writeLines('Structured Random Forest')
      varimp_rf_structured %>% arrange(desc(Overall))
      
      varimp_glm_structured = varImp(mod_glm_structured)
      writeLines('Structured GLM')
      varimp_glm_structured %>% arrange(desc(Overall))
      
      
      mod_rf_all = 
        randomForest(
          richness ~ 
            tolerance_steepness + 
            noise_level + 
            seed_output + 
            dispersal_neighborhood +
            age_structure, 
          importance = TRUE,
          data = data_richness
        )
      
      mod_glm_all = 
        glm(
          richness ~ 
            tolerance_steepness + 
            noise_level + 
            seed_output + 
            dispersal_neighborhood +
            age_structure, 
          family = gaussian,
          data = data_richness
        )
      
      varimp_rf_all = varImp(mod_rf_all)
      writeLines('All Data Random Forest')
      varimp_rf_all %>% arrange(desc(Overall))
      
      varimp_glm_all = varImp(mod_glm_all)
      writeLines('All Data GLM')
      varimp_glm_all %>% arrange(desc(Overall))
      
      mod1 = 
        glm(
          richness ~ tolerance_steepness + noise_level + mean_fecundity + dispersal_neighborhood, 
          data = data_richness
        )
      mod2 = 
        glm(
          richness ~ tolerance_steepness + noise_level + mean_fecundity + dispersal_neighborhood + age_structure, 
          data = data_richness
        )
      
      ## Testing whether the reduction in deviance by including age structure is significant using likelihood ratio test
      writeLines('Testing whether the reduction in deviance by including age structure is significant using likelihood ratio test')
      print(anova(mod1, mod2, test = 'LRT'))
      
      dtf_practice =
        tibble(
          predicted_richness = predict(mod_glm_unstructured, data = train_data),  
          observed_richness = train_data %$% richness
        )
      
      plot_prediction_practice = 
        dtf_practice %>%
        ggplot(aes(predicted_richness, observed_richness)) +
        geom_point() +
        geom_abline(aes(slope = 1, intercept = 0), color = rgb(153 / 255, 0, 0), size = 1) +
        geom_smooth(method = 'lm', color = 'blue') +
        ggtitle('A') +
        theme(aspect.ratio = 1) +
        labs(
          x = 'predicted richness',
          y = 'observed richness'
        ) +
        theme(axis.title = element_text(size = 18))
      
      writeLines('Predicted vs Observed richness on training (unstructured) data')
      summary(lm(observed_richness ~ predicted_richness, data = dtf_practice))
      
      dtf =
        tibble(
          predicted_richness = predict(mod_glm_unstructured, data = new_data),  
          observed_richness = new_data %$% richness
        )
      
      plot_prediction = 
        dtf %>%
        ggplot(aes(predicted_richness, observed_richness)) +
        geom_point() +
        geom_abline(aes(slope = 1, intercept = 0), color = rgb(153 / 255, 0, 0), size = 1) +
        geom_smooth(method = 'lm', color = 'blue') +
        ggtitle('B') +
        theme(aspect.ratio = 1) +
        labs(
          x = 'predicted richness',
          y = 'observed richness'
        ) +
        theme(axis.title = element_text(size = 18))
      
      writeLines('Predicted vs Observed richness on age-structured data')
      summary(lm(observed_richness ~ predicted_richness, data = dtf))
      
      
      gridExtra::grid.arrange(plot_prediction_practice, plot_prediction, nrow = 1)
      
    }
    
    ## Compare species abundances between structured and unstructured data
    if(TRUE){
      setwd('~/DispersalLimitation/Data/20201203_R2/')
      final_state = get(load('final_state_all_scenarios.RData'))
      
      abundances = 
        final_state %>%
        group_by(
          age_structure,
          tolerance_steepness,
          noise_level,
          mean_fecundity,
          dispersal_neighborhood,
          species
        ) %>%
        count %>%
        ungroup %>%
        mutate(n = n / 10)
      
      result_wider = 
        abundances %>%
        filter(tolerance_steepness > 0, is.finite(dispersal_neighborhood)) %>%
        mutate(age_structure = ifelse(age_structure, 'structured', 'unstructured')) %>%
        pivot_wider(names_from = age_structure, values_from = n) %>%
        replace_na(list(structured = 0, unstructured = 0))
      
      gof_test = 
        result_wider %>%
        filter(structured * unstructured > 0) %>%
        group_by(
          tolerance_steepness,
          noise_level,
          mean_fecundity,
          dispersal_neighborhood
        ) %>%
        summarize(
          pval = chisq.test(x = structured, p = unstructured, rescale.p = TRUE)$p.value,
          .groups = 'drop')
      
      
      plot_abundance_predictions = 
        result_wider %>% 
        mutate(
          noise_level = c('smooth', 'noisy')[match(noise_level, c(0, .5))], 
          noise_level = factor(noise_level, levels = c('smooth', 'noisy')),
          mean_fecundity = as.factor(mean_fecundity),
          dispersal_scale = round(dispersal_neighborhood / average_variogram_range95, 1),
          tolerance_steepness = c('tol. grad. = gradual', 'tol. grad. = steep')[match(tolerance_steepness, c(15, 100))]
        ) %>% 
        ggplot(aes(round(unstructured), round(structured), color = mean_fecundity, shape = noise_level)) + 
        geom_point() + 
        geom_abline(aes(slope = 1, intercept = 0)) + 
        facet_grid(
          tolerance_steepness ~ dispersal_scale,
          labeller = label_bquote(cols = paste('relative disp. scale = ', .(dispersal_scale)))
        ) +
        scale_x_log10() +
        scale_y_log10() +
        labs(
          x = 'abundances: unstructured',
          y = 'abundances: age-structured',
          shape = 'local \n environment',
          color = 'mean \n fecundity'
        ) +
        theme(
          strip.background = element_rect(fill = yellow),
          aspect.ratio = 1
        )
      
      result_wider %>%
        filter(structured * unstructured > 0) %>%
        group_by(
          tolerance_steepness,
          dispersal_neighborhood
        ) %>%
        summarize(
          R2 = summary(lm(structured ~ unstructured))$r.squared, 
          .groups = 'drop'
        ) %>%
        pivot_wider(names_from = dispersal_neighborhood, values_from = R2)
    }
    
  }
  
  if(Fig_niche_overlap){
    ind = 
      scenarios %>% 
      filter(
        replicate > 0 &
        dispersal_neighborhood == 1 & 
        mean_fecundity == 100 & 
        noise_level == 0 &
        age_structure == FALSE &
        neutrality == FALSE
      ) %>% 
      pull(scenario)
    
    files = paste0('scenario_', ind, '.RData')
    
    dtf = NULL
    for(char in files){
      dat = get(load(char))

      dtf =
        dtf %>%
        bind_rows(
          tibble(
            dat$scenario,
            dat$site_info,
            species = unlist(dat$species_map[, 201])
          ) %>%
            left_join(dat$fecundity, by = 'species') %>%
            left_join(dat$s50, by = 'species')
        ) %>%
        select(tolerance_steepness, replicate, site, x, y, stress, species, fecundity, s50)
    }

    dtf %<>%
      mutate(regime = ifelse(tolerance_steepness == 15, 'gradual', 'steep')) %>%
      filter(species > 0)
    
    abundance_by_stress = 
      dtf %>%
      group_by(regime, replicate, stress, species, fecundity) %>%
      count %>%
      ungroup
    
    overlap = 
      ddply(
        abundance_by_stress, 
        .(regime, replicate), 
        function(df){
          g = 
            df %>% 
            select(stress, species, n) %>%
            pivot_wider(
              names_from = species,
              values_from = n
            )
          
          g[is.na(g)] = 0
          
          h = 0
          for(i in seq(ncol(g) - 1)) for(j in (i + 1):ncol(g))
            h = c(h, sum(g[, i] * g[, j]) / sum(g[, i]) / sum(g[, j]))

          mean(h)
        }
      ) %>%
      as_tibble %>%
      rename(overlap = V1)
    
    x = overlap %>% filter(regime == 'gradual') %$% overlap
    y = overlap %>% filter(regime == 'steep') %$% overlap
    
    print(t.test(x, y))
    
    data = 
      dtf %>%
      filter(replicate == 10)
    
    foo = 
      data %>%
      group_by(regime, stress, species) %>%
      count %>%
      ungroup 
    
    bar = 
      data %>%
      group_by(regime, species) %>%
      count %>%
      ungroup %>%
      rename(abundance = n) %>%
      group_by(regime) %>% 
      mutate(rank = dense_rank(desc(abundance))) %>%
      ungroup
    
    data_for_plots =
      data %>%
      left_join(bar) %>%
      mutate(logf = log10(1e4) - log10(fecundity)) ## defining logf as the log of the seed size (not fecundity)
    
    densities = NULL
    for(thereg in c('gradual', 'steep')){
      dtf = 
        data_for_plots %>%
        filter(regime == thereg) %$%
        logf %>%
        density(bw = .08)
      
      dens = 
        tibble(
          regime = thereg, 
          logf = dtf$x, 
          density = dtf$y
        ) %>%
        mutate(
          density = 
            density / mean(density) * mean(bar %>% filter(regime == thereg) %$% abundance)
        )
      
      densities = rbind(densities, dens)
    }
    
    densities %<>% 
      group_by(regime) %>% 
      mutate(peak = splus2R::peaks(density)) %>%
      ungroup
    
    zz = NULL
    for(thereg in c('gradual', 'steep')){
      peaks_reg = 
        densities %>% 
        filter(peak == TRUE, regime == thereg)
      
      foo = 
        data_for_plots %>% 
        filter(regime == thereg) %>% 
        select(species, logf) %>% 
        unique %>% 
        arrange(logf)
      
      peak_species = 
        as.numeric(sapply(peaks_reg$logf, function(peak) with(foo, species[which.min(abs(logf - peak))])))
      
      zz = 
        zz %>%
        bind_rows(
          tibble(
            regime = thereg,
            species = peak_species
          )
        )
      
    }
    
    plot_stems = 
      data_for_plots %>%
      ggplot(aes(logf, abundance)) +
      geom_segment(aes(x = logf, y = abundance - abundance, xend = logf, yend = abundance)) +
      geom_point(color = rgb(153 / 255, 0, 0)) +
      facet_wrap(~regime) +
      theme(legend.position = 'none') +
      theme(strip.background = element_rect(fill = "#FFCD46")) +
      theme(aspect.ratio = 1) +
      xlab(expression(paste(log[10],'(fecundity)'))) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      scale_x_continuous(breaks = 1:6 / 2)
    
    
    plotA = 
      data_for_plots %>%
      ggplot(aes(logf, abundance)) +
      geom_segment(aes(x = logf, y = abundance - abundance, xend = logf, yend = abundance)) +
      geom_point(color = rgb(153 / 255, 0, 0)) +
      geom_line(aes(logf, density), data = densities, color = 'darkgrey') +
      facet_wrap(~regime) +
      theme(legend.position = 'none') +
      theme(strip.background = element_rect(fill = "#FFCD46")) +
      theme(aspect.ratio = 1) +
      xlab(expression(paste(log[10],'(seed size)'))) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) 
      #scale_x_continuous(breaks = 1:6 / 2)
    
    plotB = 
      data_for_plots %>%
      filter(regime == 'gradual') %>%
      filter(species %in% (zz %>% filter(regime == 'gradual') %$% species)) %>%
      bind_rows(
        data_for_plots %>%
          filter(regime == 'steep') %>%
          filter(species %in% (zz %>% filter(regime == 'steep') %$% species))
      ) %>%
      mutate(species = factor(species)) %>%
      ggplot(aes(stress, group = species, color = species, fill = species, alpha = 0)) +
      geom_density(bw = bw.ucv(x = seq(1,50, l = 2e3), nb = 200)) +
      # geom_density() +
      facet_wrap(~regime) +
      theme(legend.position = 'none') +
      theme(strip.background = element_rect(fill = "#FFCD46")) +
      theme(aspect.ratio = 1) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    gridExtra::grid.arrange(plotA, plotB)
    
    plot_mass = 
      data_for_plots %>%
      ggplot(aes(-logf, abundance)) +
      geom_segment(aes(x = -logf, y = abundance - abundance, xend = -logf, yend = abundance)) +
      geom_point(color = rgb(153 / 255, 0, 0)) +
      geom_line(aes(-logf, density), data = densities, color = rgb(153/255, 0, 0)) +
      facet_wrap(~regime) +
      theme(legend.position = 'none') +
      theme(strip.background = element_rect(fill = "#FFCD46")) +
      theme(aspect.ratio = 1) +
      xlab(expression(paste(log[10],'(seed size)'))) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      scale_x_continuous(breaks = 1:6 / 2)
    
  }
  
}

