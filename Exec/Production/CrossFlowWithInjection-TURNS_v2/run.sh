sh clean.sh
make -j12
mpirun -np 12 ./PeleLMeX2d.gnu.MPI.ex input.2d-regt| tee case_output.log
#python plot_pelelmx_velocity.py plt00000 -o velocity.png
#eog velocity.png_mag.png

#python plot_temperature.py plt00000 
#python plot_ch4.py plt00000 
#python plot_streamlines.py #plt00000 -o streamlines.png
#eog velocity_streamlines.png
#eog CH4.png
#eog temp.png

#python plot_time_temp_hr.py
#eog output/plt00000_temp_hr_mesh.png

