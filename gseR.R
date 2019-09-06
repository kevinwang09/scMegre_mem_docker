library(googleComputeEngineR)
project = "scmerge"
zone = "australia-southeast1-a"

gce_global_project(project)
gce_global_zone(zone)
# gce_get_project()
# gce_list_zones(project)
# View(gce_list_machinetype()$items)

(tag = "gcr.io/scmerge/scmerge_mem_docker:biocsing")

vm <- gce_vm(template = "rstudio", 
             name = "scmergebiocsing",
             disk_size_gb = 10,
             predefined_type = "n1-standard-1",
             dynamic_image = tag,
             user = "rstudio", 
             password = "pushu")


vm <- gce_ssh_setup(vm,
                    username = "rstudio",
                    key.pub = "~/.ssh/id_rsa.pub",
                    key.private = "~/.ssh/id_rsa")
