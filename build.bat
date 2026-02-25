@REM gfortran -O2 -fdefault-real-8 -Wall -ffree-line-length-none hoch.f90 common/micro_time.f90 common/linalg.f90 common/json.f90 common/jsonx.f90 common/database_m.f90 common/udp_windows_m.f90 common/connection_m.f90 vehicle.f90 sim.f90 main.f90 -o main.exe -lws2_32
gfortran -fdefault-real-8 ^
-ffree-line-length-none ^
hoch.f90 ^
common/micro_time.f90 ^
common/linalg.f90 ^
common/json.f90 ^
common/jsonx.f90 ^
common/database_m.f90 ^
common/udp_windows_m.f90 ^
common/connection_m.f90 ^
controller.f90 ^
vehicle.f90 ^
sim.f90 ^
main.f90 ^
-o main.exe -lws2_32