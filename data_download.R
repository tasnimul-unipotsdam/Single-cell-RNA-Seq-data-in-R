# Script to download the 20k NSCLC DTC dataset from 10x Genomics

# 1. Create directory if it doesn't exist
if (!dir.exists("D:/RNA-Seq/data")) {
  dir.create("D:/RNA-Seq/data")
}

# 2. Define URL and Destination
data_url <- "https://cf.10xgenomics.com/samples/cell-exp/6.1.2/20k_NSCLC_DTC_3p_nextgem_Multiplex/20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5"
dest_file <- "D:/RNA-Seq/data/20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5"

# 3. Download the file
print("Downloading dataset... this may take a few minutes.")
download.file(url = data_url, 
              destfile = dest_file, 
              mode = "wb") # 'wb' is important for binary files like .h5

print("Download complete. File saved to ../data/")