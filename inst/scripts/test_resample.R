
root <- "/mnt/storage/data/edna/mothur/projects/subsample"

filepairs = auntie::list_filepairs(file.path(root, "raw"))
default_path <- file.path(root, "test_sub")
charlier::make_path(default_path)
ok <- auntie::bbmap_subsample(filepairs, path = default_path)

fp <- auntie::list_filepairs(default_path)
ids <- auntie::fastq_ids(fp)