rm PeleLMeX3d.gnu.MPI.ex
make -j12
mpirun -np 12 ./PeleLMeX3d.gnu.MPI.ex input.3d-regt| tee case_output.log
python plot_all.py .
python plot_mesh.py .


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

