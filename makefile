sims_output = ./output/sims_out/losses_candes2.csv \
              ./output/sims_out/losses_candes.csv \
              ./output/sims_out/losses_em.csv \
              ./output/sims_out/losses_soft.csv \
              ./output/sims_out/losses_stein.csv \
              ./output/sims_out/losses_trunc.csv

mean_theta = ./output/generated_thetas/theta_candes2.txt \
             ./output/generated_thetas/theta_candes.txt \
             ./output/generated_thetas/theta_em.txt \
             ./output/generated_thetas/theta_soft.txt \
             ./output/generated_thetas/theta_stein.txt \
             ./output/generated_thetas/theta_trunc.txt

all : change_sv sims nba

## Create changing core figure in paper.
.PHONY : change_sv
change_sv : ./code/changing_core.R
	mkdir -p output/figures
	Rscript ./code/changing_core.R

## Run simulations
$(sims_out) $(mean_theta) : ./code/big_sims.R
	mkdir -p output/generated_thetas
	mkdir -p output/sims_out
	Rscript ./code/big_sims.R

## Plot simulations results
.PHONY : sims
sims : $(sims_out) $(mean_theta)
	mkdir -p output/figures
	Rscript ./code/plot_sims.R

## Run NBA Analysis
.PHONY : nba
nba : ./code/east_analyze.R ./code/west_analyze.R
	mkdir -p output/nba_results
	Rscript ./code/east_analyze.R
	Rscript ./code/west_analyze.R
