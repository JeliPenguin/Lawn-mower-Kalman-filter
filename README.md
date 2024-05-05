# COMP0130-CW1: Integrated Navigation for Robotic Lawnmower

This project aims to utilize simulated sensor data to compute the best possible horizontal position, velocity, and heading solution for a robotic lawnmower operating in a simulated environment in London. The sensors include a GNSS receiver, wheel speed sensors, a magnetic compass, and a low-cost MEMS gyroscope.

## Repository Structure


- `/COMP0130 Coursework 1 instructions 2024.pdf` - Contains the instructions for this coursework.
- `/report` - Includes the project report and documentation for the code.
- `/code` - Contains all CSV files for sensor measurements, as well as the MATLAB source code.

## Sensor Specifications

- GNSS Receiver: Assumed error specifications detailed in coursework instructions.
- Wheel Speed Sensors: Incorporates errors from varying tyre radii and wheel slip.
- Magnetic Compass and Gyroscope: Utilized for heading measurements with specified error characteristics.

## Getting Started

1. Clone the repository to your local machine.
2. Ensure you have MATLAB installed (the project is developed using MATLAB).
3. Navigate to the `/code` directory and run the main script to generate navigation solutions.