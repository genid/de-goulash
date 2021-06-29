#!/usr/bin/env python3
# -*- coding: utf-8 -*-


## clean temporal files generated during single cell pipeline

import os
import subprocess

if __name__ == "__main__":

    dirpath = os.sys.argv[1]

    target_folder = dirpath + "/iter2"
    target_file = [file_ for file_ in os.listdir(target_folder) if file_.endswith("merge_clusters.vcf")]
    if len(target_file) <= 0: # nothing to delete
        print("INFO: There are no temporal files to remove!")
        os.sys.exit()

    target = target_folder + "/" + target_file[0]
    if os.path.isfile(target):
        tmp = "{} {} {} {}".format(dirpath + "/cells", dirpath + "/data",
                               dirpath + "/iter1/pu", dirpath + "/iter2/pu")
        cmd = "rm -rf {}".format(tmp)
        subprocess.check_output(cmd, shell=True)
        tmp_dir = dirpath + "/iter1/"
        tmp = "{}*bam {}*bam.bai {}*vcf* {}*csv*".format(tmp_dir, tmp_dir, tmp_dir, tmp_dir)
        cmd = "rm {}".format(tmp)
        subprocess.check_output(cmd, shell=True)
        tmp_dir = dirpath + "/iter2/"
        tmp = "{}cells_merge_*".format(tmp_dir)
        cmd = "rm {}".format(tmp)
        subprocess.check_output(cmd, shell=True)

    target_dir = dirpath + "/output/"
    cmd = 'find {} -name "*.vcf*" -delete'.format(target_dir)
    subprocess.check_output(cmd, shell=True)
