################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_UPPER_SRCS += \
../ADTs/FEnv.C \
../ADTs/File.C \
../ADTs/Global.C 

OBJS += \
./ADTs/FEnv.o \
./ADTs/File.o \
./ADTs/Global.o 

C_UPPER_DEPS += \
./ADTs/FEnv.d \
./ADTs/File.d \
./ADTs/Global.d 


# Each subdirectory must supply rules for building sources it contributes
ADTs/%.o: ../ADTs/%.C
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -D_FILE_OFFSET_BITS=64 -I.././ADTs -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


