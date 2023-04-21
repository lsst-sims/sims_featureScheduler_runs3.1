python ../starter/starter.py --outDir ../starter
python ../filter_swap/filters_uy.py --outDir ../filter_swap
python ../filter_swap/filters_uz.py --outDir ../filter_swap
python ../noroll/noroll.py --outDir ../noroll
# Can't do these because they need ddf_grid.
# python ../roll_back_start/old_mjd0.py --outDir ../roll_back_start
# python ../roll_back_start/old_mjd0_uz.py --outDir ../roll_back_start
python ../threeseasonroll/threeroll.py --outDir ../threeseasonroll
python ../weather/weather.py --cloud_offset_year 0  --outDir ../weather
python ../weather/weather.py --cloud_offset_year 1  --outDir ../weather
python ../weather/weather.py --cloud_offset_year 2  --outDir ../weather
python ../weather/weather.py --cloud_offset_year 4  --outDir ../weather
python ../weather/weather.py --cloud_offset_year 6  --outDir ../weather
python ../weather/weather.py --cloud_offset_year 8  --outDir ../weather
python ../weather/weather.py --cloud_offset_year 10 --outDir ../weather
python ../weather/weather.py --cloud_offset_year 12 --outDir ../weather
python ../weather/weather.py --cloud_offset_year 14 --outDir ../weather
python ../weather/weather.py --cloud_offset_year 16 --outDir ../weather
python ../weather/weather.py --cloud_offset_year 18 --outDir ../weather
python ../weather/weather.py --cloud_offset_year 20 --outDir ../weather
python ../weather/weather.py --cloud_offset_year 30 --outDir ../weather
python ../weather/weather.py --cloud_offset_year 31 --outDir ../weather
python ../weather/weather.py --cloud_offset_year 35 --outDir ../weather
python ../weather/weather.py --cloud_offset_year 36 --outDir ../weather
python ../selfcal_footprints/heavy_plane.py --outDir ../selfcal_footprints
python ../selfcal_footprints/light_plane.py --outDir ../selfcal_footprints
python ../selfcal_footprints/no_plane.py    --outDir ../selfcal_footprints 
python ../upu/more_u.py --utw 12  --outDir ../upu 
python ../upu/more_u.py --utw 24  --outDir ../upu 
python ../upu/more_u.py --utw 50  --outDir ../upu 
python ../upu/more_u.py --utw 100 --outDir ../upu
python ../upu/more_u.py --utw 150 --outDir ../upu
python ../upu/more_u.py --utw 200 --outDir ../upu
python ../upu/upu.py --outDir ../upu
