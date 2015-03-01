################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../geometry/bitmap/affine.c 

OBJS += \
./geometry/bitmap/affine.o 

C_DEPS += \
./geometry/bitmap/affine.d 


# Each subdirectory must supply rules for building sources it contributes
geometry/bitmap/%.o: ../geometry/bitmap/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -D__STDC_LIMIT_MACROS -D__STDC_FORMAT_MACROS -I.././ -I.././ -O3 -Wall -c -fmessage-length=0 -static  -static-libgcc -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


