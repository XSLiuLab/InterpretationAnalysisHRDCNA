library(sigminer)
nmfall <- readRDS("/public/home/yaohz/project/HRD/data/nmfall.rds")
print('11')
sigprofiler_extract(nmfall,
                    output = "/public/home/yaohz/project/HRD/data/sig_all",
                    range = 2:30,
                    nrun = 100,
                    init_method = "random",
                    is_exome = FALSE,
                    use_conda = F
)
