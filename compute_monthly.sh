#!/bin/sh
python compute_sg_runoff_monthly.py --year="$1" --month=1 --routing="$2" --netcdf
python compute_sg_runoff_monthly.py --year="$1" --month=2 --routing="$2" --netcdf
python compute_sg_runoff_monthly.py --year="$1" --month=3 --routing="$2" --netcdf
python compute_sg_runoff_monthly.py --year="$1" --month=4 --routing="$2" --netcdf
python compute_sg_runoff_monthly.py --year="$1" --month=5 --routing="$2" --netcdf
python compute_sg_runoff_monthly.py --year="$1" --month=6 --routing="$2" --netcdf
python compute_sg_runoff_monthly.py --year="$1" --month=7 --routing="$2" --netcdf
python compute_sg_runoff_monthly.py --year="$1" --month=8 --routing="$2" --netcdf
python compute_sg_runoff_monthly.py --year="$1" --month=9 --routing="$2" --netcdf
python compute_sg_runoff_monthly.py --year="$1" --month=10 --routing="$2" --netcdf
python compute_sg_runoff_monthly.py --year="$1" --month=11 --routing="$2" --netcdf
python compute_sg_runoff_monthly.py --year="$1" --month=12 --routing="$2" --netcdf