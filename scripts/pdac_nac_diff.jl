using CSV
using ISCHIA
using DataFrames

CC = "CC07"

df_yes = DataFrame(CSV.File("outputs/pdac_nac/$(CC)_cooccurr_mat_Yes.csv"))
df_no = DataFrame(CSV.File("outputs/pdac_nac/$(CC)_cooccurr_mat_No.csv"))

sig_yes = find_differentially_cooccurring_LR_pairs(df_yes, df_no, 0.05, 0.05)
CSV.write("outputs/pdac_nac/$(CC)_differentially_cooccurring_pairs_Yes.csv", sig_yes)
println("Saved $(CC)_differentially_cooccurring_pairs_Yes.csv")

sig_no = find_differentially_cooccurring_LR_pairs(df_no, df_yes, 0.05, 0.05)
CSV.write("outputs/pdac_nac/$(CC)_differentially_cooccurring_pairs_No.csv", sig_no)
println("Saved $(CC)_differentially_cooccurring_pairs_No.csv")