# sgm-optimization

`Spectrome` is a combination of the words "spectrum" and "connectome". This package is the collection of codes that constructed the analysis various modeling papers. They are listed here:

1. ["Spectral graph model of brain oscillations"](https://onlinelibrary.wiley.com/doi/full/10.1002/hbm.24991)
2. ["Spectral graph theory of brain oscillations -- revisited and improved"](https://www.sciencedirect.com/science/article/pii/S1053811922000490)
3. ["Altered excitatory and inhibitory neuronal subpopulation parameters are distinctly associated with tau and amyloid in Alzheimer’s disease"](https://elifesciences.org/articles/77850#s3)
4. ["Stability and dynamics of a spectral graph model of brain oscillations"](https://direct.mit.edu/netn/article/7/1/48/112205/Stability-and-dynamics-of-a-spectral-graph-model). Refer to this paper for the stability analysis that is used in the optimization codes.
5. ["Structure-function models of temporal, spatial, and spectral characteristics of non-invasive whole brain functional imaging."](https://www.frontiersin.org/journals/neuroscience/articles/10.3389/fnins.2022.959557/full)
6. ["Impaired long-range excitatory time scale predicts abnormal neural oscillations and cognitive deficits in Alzheimer's disease "](https://www.researchsquare.com/article/rs-2579392/v3)

This repository is developed based on the original model's [repository](https://github.com/Raj-Lab-UCSF/spectrome).

The spectral graph model (SGM) is a brain structure-function model that simulates brain activity power spectrum given a structural connectome. The model is linear, low-dimensional, and provides an analytical relationship between the brain's structural and functional patterns.


## Abstract:
Mathematical modeling of the relationship between the functional activity and the structural wiring of the brain has largely been undertaken using non-linear and biophysically detailed mathematical models with regionally varying parameters. While this approach provides us a rich repertoire of multistable dynamics that can be displayed by the brain, it is computationally demanding. Moreover, although neuronal dynamics at the microscopic level are nonlinear and chaotic, it is unclear if such detailed nonlinear models are required to capture the emergent meso- (regional population ensemble) and macro-scale (whole brain) behavior, which is largely deterministic and reproducible across individuals. Indeed, recent modeling effort based on spectral graph theory has shown that an analytical model without regionally varying parameters can capture the empirical magnetoencephalography frequency spectra and the spatial patterns of the alpha and beta frequency bands accurately. 

In this work, we demonstrate an improved hierarchical, linearized, and analytic spectral graph theory-based model that can capture the frequency spectra obtained from magnetoencephalography recordings of resting healthy subjects. We reformulated the spectral graph theory model in line with classical neural mass models, therefore providing more biologically interpretable parameters, especially at the local scale. We demonstrated that this model performs better than the original model when comparing the spectral correlation of modeled frequency spectra and that obtained from the magnetoencephalography recordings. This model also performs equally well in predicting the spatial patterns of the empirical alpha and beta frequency bands.

## Set-up:

First clone the environment to your computer, either download this repo as a `.zip` file or `git clone https://github.com/Raj-Lab-UCSF/spectrome.git`.

Set up a [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html) if you do not have all the packages/compatible versions. The list of dependencies is listed in `environment.yml`.

Set-up environment using conda, detailed instructions can be found [here](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html). Or after cloning this repository, go to the repo by typing `cd spectrome` and then typing:
`conda env create -f environment.yml`

If conda complains about not finding packages/libraries, make sure `conda-forge` is in the list of channels being searched by `conda`.
You may add `conda-forge` to your list of channels with the command: `conda config --add channels conda-forge`.

The default name of the environment is `spectrome`, activate the environment with `source activate spectrome`, and deactivate with `source deactivate` or `conda deactivate`.

If you want to be able to run `spectrome` from anywhere, just add it's path to your PYTHONPATH. For instance, if you downloaded `spectrome` to `~/Documents/spectrome` do `export PYTHONPATH=$PYTHONPATH:~/Documents/spectrome`. You may have to restart your terminal to make sure this change takes effect.

After completing the set-up for conda environment and `spectrome` path, you may go to the `spectrome` folder and type `jupyter notebook` or `jupyter lab` in your terminal to run the Jupyter notebooks.

## Files:
 - `../spectrome/notebooks`: contains two jupyter notebooks, `run_model_example.ipynb` is the basic simulation of frequency spectrums with default parameters for the HCP template connectome, and `reproduce_MEG_modeled_spectra.ipynb` compares optimized modeled spectra with the MEG spectra.

 - `../spectrome/data`: contains intermediate data.
    - `mean80_fibercount/length.csv`: HCP template connectome and distance matrix.
    - `individual_connectomes_reordered.nc`: individual subject's connectivity matrices (N = 36).
    - `individual_psd_reordered_matlab.nc`: individual subject's MEG spectra (N = 36).

- To run optimization, refer to the scripts in the scripts folder. `sgm_fit_local.py` runs the local sgm for every region separately. `sgm_fit.py` runs the global sgm and estimates global model parameters. 

- The jobscripts have `.sh` files to launch optimization jobs either in and SGE or a SLURM cluster. You'll need to modify them according to the cluster you're working with. 
