# Assume larger_matrices is a list of larger matrices
# and smaller_matrices is a list of smaller matrices

# Loop through each smaller matrix
for (i in seq_along(smaller_matrices)) {
  smaller_mat <- smaller_matrices[[i]]
  
  # Extract the corresponding larger matrix
  larger_mat <- larger_matrices[[i]]
  
  # Append new columns from the larger matrix to the smaller matrix
  smaller_mat$CHRType <- larger_mat$CHRType
  smaller_mat$K <- larger_mat$K
  
  # If the smaller matrix is shorter than the larger matrix,
  # repeat the rows of the larger matrix until it matches the size of the smaller matrix
  if (nrow(smaller_mat) < nrow(larger_mat)) {
    repeats <- ceiling(nrow(smaller_mat) / nrow(larger_mat))
    larger_mat_repeated <- larger_mat[rep(1:nrow(larger_mat), each = repeats), ]
    larger_mat_repeated <- larger_mat_repeated[1:nrow(smaller_mat), ]
    
    # Append new columns from the repeated larger matrix to the smaller matrix
    smaller_mat$CHRType <- larger_mat_repeated$CHRType
    smaller_mat$K <- larger_mat_repeated$K
  }
  
  # Update the smaller matrix in the list
  smaller_matrices[[i]] <- smaller_mat
}