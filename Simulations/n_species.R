require(MathBioSim)
library(foreach)
library(doParallel)

linespec = function(from, to, count) {
  return((0:(count - 1)) / count * (to - from) + from)
}

file_name = "dset1"

area_length = 2
cell_count_x = 50
x_grid_points = 1000
epsilon = 1e-12

SPECIES_COUNT = 2;
B = 0.4
D = 0.2
DD = 0.001
#SIGMA_M = 0.06
LENGTH_W = 1
#MeanDensity = 50;
initial_population = 100
Intensity = 20000
Iterations = 100

N = 10

params<-list(
    "species_count"=SPECIES_COUNT,
    "area_length_x"=area_length,
    "cell_count_x"=cell_count_x,
    "spline_precision"=1e-6,
    "periodic"=FALSE
)

#numbers = sample.int(area_length*InitialDensity, SPECIES_COUNT)
#numbers = sample.int(InitialDensity, SPECIES_COUNT, replace = TRUE)

popul_x = c()
popul_s = c()

grid = linespec(0, LENGTH_W, x_grid_points)

dd_temp = c()
for (i in 1:(SPECIES_COUNT)){
    params[paste("b", i, sep = "_")] = B
    params[paste("d", i, sep = "_")] = D
    for (k in 1:(SPECIES_COUNT)){
        dd_temp = c(dd_temp, DD)
    }
    
    grid1 = linespec(0, 5*0.04, x_grid_points)
    #grid = linespec(0, length, x_grid_points)
    if(i == 1) {
        params[[paste("birth_kernel", "y", i, sep = "_")]] = dnorm(grid1, sd = 0.04)
        params[[paste("birth_kernel", "r", i, sep = "_")]] = 5*0.04
        initial_population = 100
        params[[paste("init_density", "1", sep = "_")]] = 50
    } else if(i == 2) {
        params[[paste("birth_kernel", "y", i, sep = "_")]] = dnorm(linespec(0, 5*0.06, x_grid_points), sd = 0.06)
        params[[paste("birth_kernel", "r", i, sep = "_")]] = 5*0.06
        initial_population = 200
        params[[paste("init_density", "2", sep = "_")]] = 100
    } else {
        break
    }


    




    for (j in 1:SPECIES_COUNT) {
        params[paste("dd", i, j, sep = "_")] = DD
    }
}

dd = matrix(
    dd_temp,
    nrow = SPECIES_COUNT,
    ncol = SPECIES_COUNT
)
params[["spline_precision"]] = 1e-6


dat_file_path = paste(file_name, "dat", "csv", sep = ".")
dat_file<-file(dat_file_path);
header = ""
header = paste(header, "sd_w_in",";")
header = paste(header, "sd_w_ex",";")

for (k in 1:(SPECIES_COUNT)){
    header = paste(header, "N",k,";")
}

writeLines(header, dat_file);

close(dat_file);

write_data = function(sd_w_in, sd_w_ex, N1, N2) {
            data = c()
            data = c(data, sd_w_in)
            data = c(data, sd_w_ex)
            data = c(data, N1)
            data = c(data, N2)
            cat(data, sep = ";", file = dat_file_path, append = TRUE, fill = FALSE);
            cat('\n', file = dat_file_path, append = TRUE, fill = FALSE)
        }

seed_counter = 1
for (sd_w_ex in seq(0.05, 0.21, 0.01)){
    for(sd_w_in in seq(0.2, 0.21, 0.01)) {
        print(paste("sd_w_in: ", sd_w_in, "; sd_w_ex: ", sd_w_ex))

        grid_in = linespec(0, 5*sd_w_in, x_grid_points)
        grid_ex = linespec(0, 5*sd_w_ex, x_grid_points)
        params[[paste("death_kernel", "y", 1, 1, sep = "_")]] = dnorm(grid_in, sd = sd_w_in)
        params[[paste("death_kernel", "y", 1, 2, sep = "_")]] = dnorm(grid_ex, sd = sd_w_ex)
        params[[paste("death_kernel", "y", 2, 1, sep = "_")]] = dnorm(grid_ex, sd = sd_w_ex)
        params[[paste("death_kernel", "y", 2, 2, sep = "_")]] = dnorm(grid_in, sd = sd_w_in)
        
        params[[paste("death_kernel", "r", 1, 1, sep = "_")]] = 5*sd_w_in
        params[[paste("death_kernel", "r", 1, 2, sep = "_")]] = 5*sd_w_ex
        params[[paste("death_kernel", "r", 2, 1, sep = "_")]] = 5*sd_w_ex
        params[[paste("death_kernel", "r", 2, 2, sep = "_")]] = 5*sd_w_in

        #print(params)

        time_start<-Sys.time();
        N1SUM = 0
        N2SUM = 0
        
        for (iter in 1:(N)){

            params[[paste("seed")]] = seed_counter
            seed_counter = seed_counter + 1

            sim = new(poisson_1d_n_species, params);
            for (k in 1:(Iterations)){
                sim$run_events(Intensity)
                print(sim$total_population)
            }
            N1SUM = N1SUM + sim$total_population[1]
            N2SUM = N2SUM + sim$total_population[2]
            print("---------")             
        }
        write_data(sd_w_in, sd_w_ex, N1SUM/N, N2SUM/N)
    }
}

#print(params)





#for (x in seq(0.009, 0.01, by=0.001)) {
#    params["dd_1_2"] = x
#    params["dd_2_1"] = x
#    do_simulation(paste("TEST", x, sep = "_"), params, 100, 100, FALSE)
#}

#print(paste("Start: ", file_name))
#dis_file_path = paste(file_name, "dis", "csv", sep = ".")

#total_time = 0
#file.create(dis_file_path)


#time_start<-Sys.time();


#time_finish<-Sys.time();
#total_time = total_time + time_finish - time_start;
#print(total_time);