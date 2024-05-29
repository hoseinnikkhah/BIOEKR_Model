![GitHub User's stars](https://img.shields.io/github/stars/hoseinnikkhah?style=social)
# Debugning Stages
* Run the code in MATLAB

***
## First issue

Several errors come up which were expected, the first thing we've encountered is in this line of code
```matlab
s_H(:,1) = (z_H^2)*v_H*G_H(:,1);
s_OH(:,1) = (z_OH^2)*v_OH*G_OH(:,1);
s_C(:,1) = (z_C^2)*v_C*G_C(:,1);
Sigma(:,1) = (F^2)*(s_H + s_OH + s_C);
```
As we are only changing values for all x nodes in initial time node we neede to add `(:,1)` to fix the error, ideally we will have the following code

```matlab
s_H(:,1) = (z_H^2)*v_H*G_H(:,1);
s_OH(:,1) = (z_OH^2)*v_OH*G_OH(:,1);
s_C(:,1) = (z_C^2)*v_C*G_C(:,1);
Sigma(:,1) = (F^2)*(s_H(:,1) + s_OH(:,1) + s_C(:,1));
```
## Second issue

In this line we are using `alpha` which is not defined in code and it should be 1alpha_C` instead. not to mention all concentration should be for Carbon or `G_C` in inital code I have used `G` which is not mentioned in this code.
```matlab
G_C(i,m+1) = G_C(i,m) + sub(i,m)*(alpha*(G_C(i+1,m) -2*G_C(i,m) + G_C(i-1,m)) + beta_C*(G(i+1,m) - G(i-1,m)) + R_C(i,m)/R_D);
```
So new code will look like this:
```matlab
G_C(i,m+1) = G_C(i,m) + sub(i,m)*(alpha_C*(G_C(i+1,m) -2*G_C(i,m) + G_C(i-1,m)) + beta_C*(G_C(i+1,m) - G_C(i-1,m)) + R_C(i,m)/R_D);
```
