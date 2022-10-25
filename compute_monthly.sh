#!/bin/sh
python compute_sg_runoff_monthly.py --year="$1" --month=1
python compute_sg_runoff_monthly.py --year="$1" --month=2
python compute_sg_runoff_monthly.py --year="$1" --month=3
python compute_sg_runoff_monthly.py --year="$1" --month=4
python compute_sg_runoff_monthly.py --year="$1" --month=5
python compute_sg_runoff_monthly.py --year="$1" --month=6
python compute_sg_runoff_monthly.py --year="$1" --month=7
python compute_sg_runoff_monthly.py --year="$1" --month=8
python compute_sg_runoff_monthly.py --year="$1" --month=9
python compute_sg_runoff_monthly.py --year="$1" --month=10
python compute_sg_runoff_monthly.py --year="$1" --month=11
python compute_sg_runoff_monthly.py --year="$1" --month=12