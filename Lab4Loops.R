#Lab4Loops.R Simulating breathing
#updating function f(c(t))=(1-q)*c(t)+q*y
#C(t)= concentration before one breath
#C(t+1)= conc before next breath
#q=fraction of air exchanged (W/V= volume exchanged over total volume
#y=gamma= external concentration
#at equilibrium C*=y
m=10
C=8
y=10
q=.3
update= function(C){(1-q)*C+q*y}
for (i in 1:m){
C=(1-q)*C+q*y
print(C)
}
