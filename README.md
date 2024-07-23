# 5G power optimization
MATLAB implementation of the power optimization in 5G networks with Massive MIMO technique using the Dinkelbach algorithm and Water Filling, both uplink and downlink, on Sum Rate and Energy Efficiency.

<br>
<img src="./images/Users.png"> <br>

# File organization

The repository is divided into two different folders. In order to start the simulation open each single folder in MATLAB and run the corresponding starting point:

- **uplink**: optimization in uplink comunication. The starting point is the file *main_uplink.m*. The problem is divided into:
    - **Sum Rate optimization**: Each user can transmit to its maximum power. <br>
    <img src="./images/SumRate-uplink.png"> <br>
    - **Energy Efficiency optimization**: the Dinkelback algorithm optimaze the ratio of the Energy Efficiency.<br>
    <img src="./images/EnergyEfficiency-uplink.png"> <br>
    <br>
    The algorithm maxime the auxilary function of the fractional problem:
    <br>
    <img src="./images/AuxiliaryFunction.png"> <br>
    The algorithm:
    <br>
    <img src="./images/Dinkelback.png"> <br>

- **downlink**: optimization in downlink comunication. The starting point is the file *main_downlink.m*. The problem is divided into:
    - **Sum Rate optimization**: The radio base station have to transmit towards every user, considering that the sum of all the power allocated to the users have to sum up the maximum power available to the station. <br>
    <img src="./images/SumRate-downlink.png"> <br>
    The power is allocated with the Water Filling algorithm:
    <br>
    <img src="./images/WaterFilling.png"> <br>
    - **Energy Efficiency optimization**: as for the uplink problem, the Dinkelback algorithm optimaze the ratio of the Energy Efficiency.<br>
    <img src="./images/EnergyEfficiency-downlink.png"> <br>

# Simulation
Each problem compare the optimal power allocation with the uniform one, considering also the case with and without interference.

The channel parameters are extracted from a normal distribution.

**uplink**:
<br>
<img src="./images/SumRate-power-uplink.png"> <br>
<img src="./images/EnergyEfficiency-power-uplink.png"> <br>

**downlink**:
<br>
<img src="./images/SumRate-power-downlink.png"> <br>
<img src="./images/WaterFilling-allocation.png"> <br>
<img src="./images/EnergyEfficiency-power-downlink.png"> <br>

