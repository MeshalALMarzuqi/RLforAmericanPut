# RLforAmericanPut
Reinforcement Learning for American Put Options
Pricing derivatives with American features remains a challenging and an interesting problem
for the financial engineering community, where standard numerical methods such as finite
differences become intractable in a multi-factor environment. In this thesis we compare the
performance of three state-of-the-art simulation based continuous Reinforcement Learning
algorithms Longstaff-Schwartz (LSM), Fitted-Q Iteration (FQI), and Least Squares Policy
Iteration (LSPI) as applied to the pricing and trading of American put options. Initially we
closely follow the paper of Schuurmans et. al. and attempt to re-create their results by
fully mathematically deriving each algorithm in order to gain in-depth understanding of they
work and by designing the necessary coding architecture, which can be built upon. The authors
report the largest payoffs given by the FQI algorithm, followed by LSPI, which gives similar
but slightly lower payoffs, followed by LSM, which is reported to give poor payoffs (sometimes
lower than its European counterpart) and very choppy exercise boundaries. Upon the initial
re-creation of the results, our payoffs agree with those of the original paper, where FQI finds
the highest payoffs, followed by LSPI and LSM. Next we seek to improve the reported results
and investigate how additional basis functions (which depend on the option’s volatility) affect
learning of the exercise policies. We find that the weights learnt for volatility basis are close
to zero, indicating that these basis were not found to be important for the continuation value
function approximation. Lastly we attempted to improve the policies learnt using the LSM
and their associated exercise boundaries, since this algorithm is the default method of choice
within the Financial Engineering community for pricing derivatives with American features. We
discover that our version of the Longstaff-Schwartz algorithm is able to gain considerably larger
payoffs and produce boundaries which are the best out of all the tested algorithms (smooth and
monotonically increasing) and we argue that this is an original contribution to the Schuurmans
et. al. paper.
