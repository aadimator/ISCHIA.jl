using CSV
using ISCHIA
using DataFrames

df1 = DataFrame(CSV.File("outputs/CC1_cooccurr_mat.csv"))
df3 = DataFrame(CSV.File("outputs/CC3_cooccurr_mat.csv"))

diff_13 = find_differentially_cooccurring_LR_pairs(df1, df3, 0.05, 0.1)

CSV.write("outputs/diff_cc1_cc3.csv", diff_13)
println("Saved to outputs/diff_cc1_cc3.csv")