################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/aux_functions.cpp \
../src/elements.cpp \
../src/export.cpp \
../src/insert_from_file.cpp \
../src/sequential_interface.cpp \
../src/smeagol.cpp \
../src/solver.cpp 

OBJS += \
./src/aux_functions.o \
./src/elements.o \
./src/export.o \
./src/insert_from_file.o \
./src/sequential_interface.o \
./src/smeagol.o \
./src/solver.o 

CPP_DEPS += \
./src/aux_functions.d \
./src/elements.d \
./src/export.d \
./src/insert_from_file.d \
./src/sequential_interface.d \
./src/smeagol.d \
./src/solver.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


