library(googleComputeEngineR)
project = "scmerge"
zone = "us-central1-a"

gce_global_project(project)
gce_global_zone(zone)
# gce_get_project()
# gce_list_zones(project)
# View(gce_list_machinetype()$items)

(tag = "gcr.io/scmerge/scmerge_mem_docker:gcbuild")

vm <- gce_vm(template = "rstudio", 
             name = "biocsing8",
             disk_size_gb = 20,
             predefined_type = "n1-standard-8",
             dynamic_image = tag,
             user = "rstudio", 
             password = "pushu")


vm <- gce_ssh_setup(vm,
                    username = "rstudio",
                    key.pub = "~/.ssh/id_rsa.pub",
                    key.private = "~/.ssh/id_rsa")
