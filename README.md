# ShapeSegmentationPOP
This project attempts to use polynomial optimization for finding a global solution to shape segmentation

The main script to be called is run_iterate.m which for the chosen number of times calls multiple_weights_fixing. This function in turn analyzes the input image and calls either write_Gams_levelset2 (for 4th degree polynomials) or write_Gams_levelset2_arplane.m (for 6th degree polynomials). These scripts in turn will create the appropriate vectors of monomials and coefficients in Gams format and call sparsePOP function (SparsePOP package needs to be set up and added to project path). SparsePOP runs an SDP solver (refer to the package's setting) and from it's returning values we'll get the coefficients of the curve that was fit.
