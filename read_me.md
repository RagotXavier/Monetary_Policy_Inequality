---
author: François Le Grand, Alaïs Martin-Baillon, Xavier Ragot
title: "Replicating package*: Should monetary policy care about
  redistribution? Optimal monetary and fiscal policy with heterogeneous
  agents.* "
date: June, 27 2024
colorlinks: true
---

This set of files replicates the graphs and data of *Should monetary
policy care about redistribution? Optimal monetary and fiscal policy
with heterogeneous agents.* Be careful of not mixing up these files with
other replicating files and they should be located in an independent
folder. Requires Matlab / Julia / Dynare. Tested on Matlab
2018b and Julia release v1.10.0.

1.  All the `.ipynb` files are [`Julia`](https://julialang.org/) notebooks.

2.  All the `.m` files can be run in Matlab or [`Octave`](https://octave.org/). Most of them also require [`Dynare`](https://www.dynare.org/).

# Quantitative assessment of the sticky-price model

## How to run?

1. Start with the computation of the steady state allocation:

*  Open `Sticky_Prices/steady_state/Main_SP.ipynb` and execute all
    cells. 
    
2. Then, simulate the dynamics of the model (IRFs and second-order moments) in one of the following cases. They are independent and can be run in any order. 

    * **Baseline calibration and the uniform truncation**  

        a. Run `Sticky_Prices/dynamics/main_SP.m` to simulate the dynamics of Economies 1 and 2.
        
        b. Run `Sticky_Prices/dynamics/main_taylor_SP.m` to simulate the dynamics of Economy 3 with the Taylor rule.
        
        c. Run `Sticky_Prices/dynamics/Do_IRFs_SP_baseline.m` to replicate Figure 1.
        
        d. Use `Sticky_Prices/steady_state/Tables_SP.ipynb` to replicate Table 7.

    *  **Baseline calibration and the refined truncation** 
        a. Open  `Sticky_Prices/dynamics/main_SP.m` and comment the line 3 `calib = 'baseline'` and uncomment the line 4 `calib = 'refined'`. 
        
        b. Run `Sticky_Prices/dynamics/main_SP.m`  to simulate the dynamics. 
        
        c. Use `Sticky_Prices/steady_state/Tables_SP.ipynb` to replicate Table 8.
        

    * **Unequal profit distribution**

        a. Run `Sticky_Prices/dynamics/unequal.m`  to simulate the dynamics. 
        
        b. Run `Sticky_Prices/dynamics/Do_IRFs_SP_unequal.m` to replicate Figure 7. 
    
3. To compute the welfare per productivity level:

    a. Run `Sticky_Prices/dynamics/welfare.m`  to simulate the dynamics of Economies 1 and 2 over 10,000 periods with different seeds and to compute the welfare differences between the two economies per level of productivity. 
    
    b. Run `Sticky_Prices/dynamics/Do_Welfare_SP.m` to replicate Figure 2.

4. To check the robustness relative to the slope of the Phillips curve:

    a. Run `Sticky_Prices/dynamics/robustness_slope_SP.m` to simulate the dynamics of the model for different values of the slope of the PC. 
    
    b. Run `Sticky_Prices/dynamics/Do_IRFs_SP_rob.m` to replicate Figure 5. 




## The details

The `Julia` file takes care of computing the steady state, while     the `.m` files simulate the dynamics of the model in the presence of aggregate shocks.

The outputs of the `Julia` file are the following `.mat` files: 

1. `todynare_SP_baseline.mat`; 
2. `todynare_SP_refined.mat`; 
3. `To_IRFs_SP_unequal.mat`. 
4. `todynare_SP_rob_20.mat` , `todynare_SP_rob_58.mat`  , `todynare_SP_rob_100.mat` , `todynare_SP_rob_263.mat`  , `todynare_SP_rob_806.mat`.


The first one is used by `main_SP.m`, `main_taylor_SP.m` and  `welfare.m` ; the second by `main_SP.m` (with the proper modification of the input file discussed above); third last one by `unequal.m`, the forth ones by `robustness_slope_SP.m`

The outputs of the previous `.m` files are respectively:

1. `To_IRFs_SP_baseline.mat`;
2. `To_IRFs_SP_taylor.mat`;
3. `To_IRFs_SP_unequal.mat`.
4. `To_IRFs_SP_rob.mat`.
5. `For_Welfare_SP.mat`.


Be careful: Depending on the input (`calib = 'baseline'` or `calib = 'refined'`), `main_SP.m` generates two different inputs that are both named `To_IRFs_SP_baseline.mat`.

The graphs can then be plotted as follows:

*  `Do_IRFs_SP_baseline.m` generates  `IRFs_SP_Eco_1_2_taylor.png` corresponding to Figure 1.
*  `Do_IRFs_SP_unequal.m` generates `IRFs_SP_uneq_Eco_1_2.png` corresponding to Figure 7. 
*  `Do_IRFs_SP_rob` generates `Robustness_Slope_SP.png`  corresponding to Figure 5.
*  `Do_Welfare_SP` generates `Welfare_SP_baseline_Eco_1_2`  corresponding to Figure 2.

Finally, the Julia file `Tables_SP.ipynb` displays Tables 7 and 8 for the first- and second-order moments of the key model variables. The file uses other outputs of the previous `.m` file, which are: `moments_eco1/2/3_baseline.mat`, `moments_taylor.mat`, `moments_eco1/2_refined.mat`.

## The steady state computation

The steady state is computed thanks to nine `Julia` notebooks. Each of
these files are commented and self-explained. Briefly:  

* `Main_SP.ipynb`: Solves the steady-state model and returns the
    truncated model (as `steady_state_dynare.mat` for `Dynare`, saved in
    the current folder). It relies on most of the other notebooks. 

* `Structures_SP.ipynb`: Structures and parameter calibration from targets.

* `Utils_SP.ipynb`: Contains some useful functions.

* `SolveAiyagari_SP.ipynb`: Solves the Aiyagari model.

* `Projection_SP.ipynb`: Computes the steady-state truncated model.

* `Projection_SP_ref.ipynb`: Computes the steady-state refined truncated
    model.

* `Ramsey_SP.ipynb`: Computes the steady-state Lagrange multipliers.

* `Simulation_SP.ipynb`: : Contains a function used to display tables of
    first and second order moments of key variables.

* `Tables_SP.ipynb`: Displays tables of first- and second-order moments of     key variables.

## Simulating the model with aggregate shocks

* The file `main_SP.m` simulates the model for the first two economies of
    the paper: Economy 1 corresponds to optimal inflation and Economy 2 to constant inflation. For the baseline truncation or the refined truncation
    depending on the choice of the input file (lines 3 and 4 telling whether `calib = 'baseline'` or `calib = 'refined'`).

*  The file `main_taylor_SP.m` simulates Economy 3 with the  Taylor rule. 
*  The file `robustness_slope_SP.m` simulates economies with different slopes of the Phillips Curve. 


* The `.m` files actually write different `Dynare` executable codes (as `.mod` files). The `.mod` files are created within the current folder and then automatically run. They are:

    * `code_dynare_baseline1.mod`, `code_dynare_baseline2.mod`, `code_dynare_taylor.mod`, which correspond to the three economies of the baseline specification;
    * `code_dynare_refined1.mod`, `code_dynare_refined2.mod`, which correspond to the two economies of the refined specification;
    * `code_dynare_unequal1.mod`, `code_dynare_unequal2.mod`, which correspond to the two economies of the specification with unequal profit distribution.
    * `code_dynare_rob1.mod`, `code_dynare_rob2.mod`, `code_dynare_ rob3.mod`, `code_dynare_rob4.mod`, `code_dynare_rob5.mod` which correspond to the economies wit the different slopes of the Phillips curve. 
    * `code_dynare_SP_welfare1.mod`, `code_dynare_SP_welfare2.mod`,  which correspond to the two economies of the baseline specification.

# Comparisons with the Reiter method

## How to run?

In this order:

1.  Open `Reiter/steady_state/Main_Reiter_Comp.ipynb` and execute all
    cells. This computes the steady-state (common to both methods)

2.  Run `Reiter/dynamics/reiter_comp.m`, which simulates the dynamics of the model with the Reiter method.

3.  Run `Reiter/dynamics/refined_comp.m`, which simulates the dynamics of the model with the refined truncation method.

4.  Run `Reiter/dynamics/code_difference.m`, which computes the differences between the outputs of the two methods. 

To compute the outputs:

* Run `Reiter/dynamics/fig_Comp_Reiter_Trunc.m` to plot  Figure 8.

* Run `Reiter/steady_state/Tables_Reiter_Comp.ipynb` (cells 1 to 4) to replicate  Table 11.

* Run `Reiter/steady_state/Tables_Reiter_Comp.ipynb` (cell 5) to replicate Table 12

## The details

The outputs of the `Julia` files computing the steady state are: `todynare_Comp_ Reiter.mat` and `todynare_Comp_Refined.mat`. These two files are used for the simulation of the dynamics of the model with the two approaches and in  `reiter_comp.m` and `refined_comp.m`, respectively. 

The outputs of the `.m` files are the following `.mat` files.

* `todiff_Reiter IRFs.mat` saves the IRFs for the Reiter method;

* `todiff_Comp_Refined.mat` saves the IRFs for the refined truncation method; 

* `moments_Comp_Reiter.mat` saves the moments for the Reiter method; 

* `moments_Comp_Refined.mat` saves the moments for the refined truncation method;

* `to_code_difference_reiter.mat` save the simulation results for the Reiter method;

* `to_code_difference_refined.mat` save the simulation results for the refined truncation method;

Those files are then processed by `code_difference.m` to generate `diff_Reiter.mat` that contains the results for comparing the two methods.  Finally:

* `Reiter/steady_dynamics/fig_Comp_Reiter_Trunc.m` plots Figure 8 called `Comp_Reiter_Trunc.png`;

* `Reiter/steady_state/Tables_Reiter_Comp.m` displays Tables 11 and 12.

# Quantitative assessment of the sticky-wage model

The procedure is very close to the one for the sticky-price model. 

## How to run?

In this order:

1. For the baseline calibration:

    a. Open `Sticky_Wages/steady_state/Main_SW.ipynb` and execute all cells. This computes the steady state. 

    b. Run `Sticky_Wages/dynamics/main_SW.m`. This simulates the dynamics of Economies 1 and 2.

    c. Run `Sticky_Wages/dynamics/main_taylor_SW.m`. This simulates the dynamics of Economy 3.

    d. Run `Sticky_Wages/dynamics/Do_IRFs_SW_baseline.m` to plot the IRFs of Figure 3. 


2. For the robustness relative to the slope of the Phillips curve:

    1. Run `Sticky_Wages/dynamics/robustness_slope_SW.m` to simulate the dynamics and replicate Figure 12. 
    
3. To compute the welfare per productivity level:

    1. Run `Sticky_Wages/dynamics/welfare_SW.m`  to simulate the dynamics of Economies 1 and 2 over 10,000 periods with different seeds and to compute the welfare differences between the two economies per level of productivity. 

    2. Run `Sticky_Wages/dynamics/Do_Welfare_SW.m` to replicate Figure 4.


## The details

They are very similar to those of the sticky-price economy. Files contains `SW` instead of `SP`. The output of the  `Julia` files computing the steady state is a file `todynare_SW.mat` that is then processed by `Dynare`. The outputs of the `Dynare`files are `To_IRFs_SW.mat` and `To_IRFs_SP_taylor.mat`, which are used by `Do_IRFs_SP_baseline.m` to plot Figure 3 called `IRFs_SW_Eco_1_4_taylor.png`.

 

## The steady state computation

The steady state computation is done via `Julia` notebooks that are similar to those of the sticky-price economy. They are suffixed by `SW` instead of `SP`.


## Simulating the model with aggregate shocks,

Again, this is very similar to the sticky-price economy.  The file `main_SW.m` simulates the model for Economy 1 (optimal inflation) and Economy 2 (constant wage inflation). The file `main_SW_taylor.m` simulates the model for Economy 3 with the Taylor Rule.
