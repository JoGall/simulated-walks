# Simulate locomotion using randomly generated X,Y-trajectories

### 1. Simulate animal locomotion using a random walk rule

Developed to model insect locomotion using a [Lévy random walk](https://en.wikipedia.org/wiki/L%C3%A9vy_flight), in conjunction with the tracking software [Ubitrail](https://github.com/JoGall/ubitrail) and analysis package [RUbitrail](https://github.com/JoGall/rubitrail). Scripts to make similar random walks found [here](https://github.com/JoGall/simulated-walks/blob/master/makeWalks.R), and several variants of random walk functions found [here](https://github.com/JoGall/simulated-walks/blob/master/walkFuns.R)

<img src="https://cloud.githubusercontent.com/assets/17113779/14388979/792ab794-fda8-11e5-98a5-a32ee04165b3.png" width="300"><img src="https://cloud.githubusercontent.com/assets/17113779/14388983/79370ec2-fda8-11e5-8ffa-f27f9e5fd525.png" width="300">

**Figure 1.** Sample 12000 frame trajectories, showing endogenous locomotion in the mealworm beetle, *Tenebrio molitor* (left), and a simulated trajectory produced by a correlated walk rule (right).

<img src="https://cloud.githubusercontent.com/assets/17113779/14388978/79282d08-fda8-11e5-9376-13fec9f40a62.png" width="300"><img src="https://cloud.githubusercontent.com/assets/17113779/14388980/792da5da-fda8-11e5-8326-58dbf3a76be9.png" width="300">
<img src="https://cloud.githubusercontent.com/assets/17113779/14388981/7930995c-fda8-11e5-85e1-9a706107da17.png" width="300"><img src="https://cloud.githubusercontent.com/assets/17113779/14388982/793347ce-fda8-11e5-9f4e-0721bad55ac8.png" width="300">

**Figure 2.** Correlated walks produced using various turning parameters, where rho = 0.1, 0.5, 0.9 and 0.994, from top left to bottom right.


### 2. Simulate movement of football players using a random walk rule

Simulates X,Y-trajectories of football players using a Lévy random random walk, with realistic velocity, acceleration and pauses. Parameters fitted in [this script](https://github.com/JoGall/simulated-walks/blob/master/fitFootballRun.R) using real player positions from a sample dataset used in [this publication](http://home.ifi.uio.no/paalh/publications/files/mmsys2014-dataset.pdf).

<img src="https://user-images.githubusercontent.com/17113779/30368517-48f38bde-9869-11e7-8cc2-4ce888136e13.png" width="400"><img src="https://user-images.githubusercontent.com/17113779/30368520-4aaed38e-9869-11e7-80c5-2c319c16754b.png" width="400">

**Figure 3.** Left: Movement of one Tromsø IL player in the first minute (1200 frames) of the [match vs. Strømsgodset](http://home.ifi.uio.no/paalh/dataset/alfheim/) (2013-11-03). Right: Simulated X,Y-trajectories for 1200 frames based upon parameters fitted to all players in this match. Background pitch created in gglot using [drawPitch.R](https://github.com/JoGall/football-xy/blob/master/drawPitch.R).
