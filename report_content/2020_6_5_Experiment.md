# Main experimental contents and conclusions of this week (China time: 2020.6.5, Friday) :

## 1.When all DCs' privacy level set to 1

### The experimental setup:

 Model: NewPregel

 DCs' privacy level: {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}

 Test graphs: Power-law graph and NetWork graph

 Sampling probability formula： P=|1/last_change| and P=|1/(rank-0.15)|

 Application: PageRank and SSSP
   
### The experimental result:
##### Under the power-law graph:
 PageRank: the average relative of all vertices error is about 22%, the number of iterations is 20.

 SSSP: the average relative of all vertices error is about 25%, the number of iterations is 57.


##### Under the network test graph:
 PageRank: the average relative of all vertices error is about 9%, the number of iterations is 20.

 SSSP: the average relative of all vertices error is about 28%, the number of iterations is 33.

## 2. Relative error among data centers,PageRank and SSSP
### The experimental setup:

 Model: NewPregel

 DCs' privacy level: {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1} and {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,20}

 Test graph: Power-law graph

 Sampling probability formula： P=|1/last_change|

 Application: PageRank and SSSP

 budget: 1

### The experimental result:
 PageRank: No matter which combination of privacy levels, no matter which test graph, the average relative error in each DC is approximately equal to the average relative error of all vertices.

 SSSP: No matter which combination of privacy levels, for Power-law Graph, the average relative error in each DC is approximately equal to the average relative error of all vertices, while for Network Graph, the average relative error is slightly higher than other DCs in DC2 and DC14.

## 3. Impact of the number of iterations
### The experimental setup:

 Model: NewPregel

 DCs' privacy level: {1111,2222,3333,4444,5555}

 Test graph: Power-law graph

 Sampling probability formula： P=|1/last_change| and P=|1/(rank-0.15)|

 Application: PageRank

 budget: 1

### The experimental result:
 Regardless of the sampling probability formula, the average relative error decreases with the increase of iteration number before the iteration number is 30.After 30 iterations, the mean relative error increases with the number of iterations. 

## 4. More uniform privacy level distribution

### The experimental setup:

 Model: NewPregel

 DCs' privacy level: {1111,2222,3333,4444,5555}

 Test graph: Power-law graph

 Sampling probability formula： P=|1/last_change| and P=|1/(rank-0.15)|

 Application: PageRank and SSSP

 budget: 1,2,3,5


### The experimental result:
#### PageRank, P=|1/last_change|: 
 Average relative error: 194%(budget=1)->117%(budget=2)
->75%(budget=3)->43%(budget=5)
 
 Wan usage: 5.9MB(budget=1)->5.66MB(budget=2)
->6.16MB(budget=3)->7.06MB(budget=5)

#### PageRank, P=|1/(rank-0.15)|: 
 Average relative error: 124%(budget=1)->111%(budget=2)
->100%(budget=3)->30%(budget=5)
 
 Wan usage: 7.32MB(budget=1)->7.06MB(budget=2)
->6.42MB(budget=3)->8.17MB(budget=5)

#### SSSP, P=|1/last_change|: 
 Average relative error: 93%(budget=1)->38%(budget=2)
->30%(budget=3)->28%(budget=5)
 
 Wan usage: 7.1MB(budget=1)->7.3MB(budget=2)
->7.3MB(budget=3)->7.2MB(budget=5)

# The next week's work
1. In view of the {1, 1,..., 1} vs. {1, 1,...1,20} these two kinds of setting. Check the each iteration's sampling rate.

2. For SSSP algorithm, when the privacy level of all DCs are set to 1, find the reason why average relative error of DC2 and DC14 are relatively large.


3. Set the budget as 0.1, that is, continue to increase the noise and observe the relative error changing with the number of iterations.
4. For PageRank, under the condition of P=|1/last_change|, the average relative error is not obvious with the change of iteration number, so the iteration number is continued to increase(>60) to observe the change of average relative error with iteration number.