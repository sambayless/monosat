################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../aiger/aiger.c 

OBJS += \
./aiger/aiger.o 

C_DEPS += \
./aiger/aiger.d 


# Each subdirectory must supply rules for building sources it contributes
aiger/%.o: ../aiger/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -DNDEBUG -I.././ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


