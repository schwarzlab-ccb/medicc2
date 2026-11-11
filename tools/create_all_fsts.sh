#!/bin/bash

python ../tools/medicc-create-cn-fst.py ../medicc/objects --max-cn 8 --sep X --wgd --prefix wgd
python ../tools/medicc-create-cn-fst.py ../medicc/objects --max-cn 8 --sep X --wgd --wgd_x2 --prefix wgd_x2
python ../tools/medicc-create-cn-fst.py ../medicc/objects --max-cn 8 --sep X --wgd --wgd_x2 --max-num-wgds 1 --prefix wgd_x2_1
python ../tools/medicc-create-cn-fst.py ../medicc/objects --max-cn 8 --sep X --wgd --total_cn --prefix wgd_total_cn
python ../tools/medicc-create-cn-fst.py ../medicc/objects --max-cn 8 --sep X --wgd --total_cn --max-num-wgds 1 --prefix wgd_total_cn_1
python ../tools/medicc-create-cn-fst.py ../medicc/objects --max-cn 8 --sep X --prefix no_wgd
python ../tools/medicc-create-cn-fst.py ../medicc/objects --max-cn 8 --sep X --wgd --max-num-wgds 1 --prefix wgd_1
python ../tools/medicc-create-cn-fst.py ../medicc/objects --max-cn 8 --sep X --wgd --max-num-wgds 2 --prefix wgd_2
python ../tools/medicc-create-cn-fst.py ../medicc/objects --max-cn 8 --sep X --wgd --prefix legacy_loh --legacy-loh
