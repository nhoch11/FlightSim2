@REM gfortran -O2 -fdefault-real-8 -Wall -ffree-line-length-none hoch.f90 common/micro_time.f90 common/linalg.f90 common/json.f90 common/jsonx.f90 common/database_m.f90 common/udp_windows_m.f90 common/connection_m.f90 vehicle.f90 sim.f90 main.f90 -o main.exe -lws2_32
gfortran -fdefault-real-8 ^
-finit-real=nan -fcheck=all -O0 -g ^
-ffree-line-length-none ^
src/hoch.f90 ^
common/micro_time.f90 ^
common/linalg.f90 ^
common/json.f90 ^
common/jsonx.f90 ^
common/database_m.f90 ^
common/udp_windows_m.f90 ^
common/connection_m.f90 ^
src/atmosphere.f90 ^
src/controller.f90 ^
src/propulsion.f90 ^
src/vehicle.f90 ^
src/sim.f90 ^
src/main.f90 ^
-o main.exe -lws2_32