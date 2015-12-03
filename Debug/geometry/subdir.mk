################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../geometry/ConvexHullCollisionDetection.cpp \
../geometry/ConvexHullDetector.cpp \
../geometry/ConvexPolygon.cpp \
../geometry/DelaunayPolypartition.cpp \
../geometry/GridHeightmap.cpp \
../geometry/PointSet.cpp \
../geometry/QuickConvexHull.cpp 

OBJS += \
./geometry/ConvexHullCollisionDetection.o \
./geometry/ConvexHullDetector.o \
./geometry/ConvexPolygon.o \
./geometry/DelaunayPolypartition.o \
./geometry/GridHeightmap.o \
./geometry/PointSet.o \
./geometry/QuickConvexHull.o 

CPP_DEPS += \
./geometry/ConvexHullCollisionDetection.d \
./geometry/ConvexHullDetector.d \
./geometry/ConvexPolygon.d \
./geometry/DelaunayPolypartition.d \
./geometry/GridHeightmap.d \
./geometry/PointSet.d \
./geometry/QuickConvexHull.d 


# Each subdirectory must supply rules for building sources it contributes
geometry/%.o: ../geometry/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++1y -D__STDC_LIMIT_MACROS -D__STDC_FORMAT_MACROS -I.././ -O0 -g3 -Wall -c -fmessage-length=0  -Wno-unused-variable -Wno-unused-but-set-variable   -Wno-sign-compare -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


