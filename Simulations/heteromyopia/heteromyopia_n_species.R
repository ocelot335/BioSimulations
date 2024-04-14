require(MathBioSim)
library(foreach)
library(doParallel)

linespec = function(from, to, count) {
  return((0:(count - 1)) / count * (to - from) + from)
}

file_name = "dset1_n_species"

area_length = 10
cell_count_x = 50
x_grid_points = 1000
epsilon = 1e-12

SPECIES_COUNT = 2;
B = 0.4
D = 0.2
DD = 0.001
Intensity = 20000
Iterations = 100
SIGMA_M1 = 0.04
SIGMA_M2 = 0.06
N01 = 100
N02 = 200

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
    
    if(i == 1) {
        params[[paste("birth_kernel", "y", i, sep = "_")]] = dnorm(linespec(0, 5*SIGMA_M1, x_grid_points), sd = SIGMA_M1)
        params[[paste("birth_kernel", "r", i, sep = "_")]] = 5*SIGMA_M1
        params[[paste("init_density", "1", sep = "_")]] = N01/area_length
    } else if(i == 2) {
        params[[paste("birth_kernel", "y", i, sep = "_")]] = dnorm(linespec(0, 5*SIGMA_M2, x_grid_points), sd = SIGMA_M2)
        params[[paste("birth_kernel", "r", i, sep = "_")]] = 5*SIGMA_M2
        params[[paste("init_density", "2", sep = "_")]] = N02/area_length
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

header = paste(header, "time1",";")
header = paste(header, "time2",";")

writeLines(header, dat_file);

close(dat_file);

write_data = function(sd_w_in, sd_w_ex, N1, N2, TIME1, TIME2) {
            pop = sim$total_population;
            data = c()
            data = c(data, sd_w_in)
            data = c(data, sd_w_ex)
            data = c(data, N1)
            data = c(data, N2)
            data = c(data, TIME1)
            data = c(data, TIME2)
            cat(data, sep = ";", file = dat_file_path, append = TRUE, fill = FALSE);
            cat('\n', file = dat_file_path, append = TRUE, fill = FALSE)
        }

seed_counter = 1
for (sd_w_ex in seq(0.01, 0.2, 0.05)){
    for(sd_w_in in seq(0.01, 0.2, 0.05)) {
        print("---")
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
        TIME1SUM = 0
        TIME2SUM = 0
        for (iter in 1:(N)){

            params[[paste("seed")]] = seed_counter
            seed_counter = seed_counter + 1

            sim = new(poisson_1d_n_species, params);
            DIED1 = 0
            DIED2 = 0
            for (k in 1:(Iterations)){
                sim$run_events(Intensity)
                popul = sim$total_population
                print(popul)
                if(popul[1] == 0 && DIED1 == 0) {
                    DIED1 = 1
                    TIME1SUM = TIME1SUM + k*Intensity
                }
                if(popul[2] == 0 && DIED2 == 0) {
                    DIED2 = 1
                    TIME2SUM = TIME2SUM + k*Intensity
                }
            }
            if(DIED1 == 0) {
                DIED1 = 1
                TIME1SUM = TIME1SUM + Iterations*Intensity*100
                }
            if( DIED2 == 0) {
                DIED2 = 1
                TIME2SUM = TIME2SUM + Iterations*Intensity*100
            }
            N1SUM = N1SUM + sim$total_population[1]
            N2SUM = N2SUM + sim$total_population[2]             
            print(paste("sd_w_in: ", sd_w_in, "; sd_w_ex: ", sd_w_ex))
        }
        write_data(sd_w_in, sd_w_ex, N1SUM/N, N2SUM/N, TIME1SUM/N, TIME2SUM/N)
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