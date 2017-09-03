function count = PoissonDist(lambda,k)
count = exp(-lambda)*lambda^k/factorial(k);