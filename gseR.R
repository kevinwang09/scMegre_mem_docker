library(googleComputeEngineR)
project = "scmerge"
zone = "australia-southeast1-a"
# zone = "asia-east2-a" ## Hong Kong server

gce_global_project(project)
gce_global_zone(zone)
# gce_get_project()
# gce_list_zones(project)
# View(gce_list_machinetype()$items)

(tag = gce_tag_container(""))

vm <- gce_vm(template = "rstudio", 
             name = "scmerge_mem",
             disk_size_gb = 10,
             predefined_type = "n1-standard-4",
             dynamic_image = tag,
             user = "rstudio", 
             password = "pushu")


vm <- gce_ssh_setup(vm,
                    username = "rstudio",
                    key.pub = "~/.ssh/id_rsa.pub",
                    key.private = "~/.ssh/id_rsa")
