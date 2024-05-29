![GitHub User's stars](https://img.shields.io/github/stars/hoseinnikkhah?style=social)
### Debugning Stages
* Run the code in MATLAB

*** First issue

Several errors come up which were expected, the first thing we've encountered is in this line of code
```
s_H(:,1) = (z_H^2)*v_H*G_H(:,1);
s_OH(:,1) = (z_OH^2)*v_OH*G_OH(:,1);
s_C(:,1) = (z_C^2)*v_C*G_C(:,1);
Sigma(:,1) = (F^2)*(s_H + s_OH + s_C);
```
As we are only changing values for all x nodes in initial time node we neede to add `(:,1)` to fix the error, ideally we will have the following code

```
s_H(:,1) = (z_H^2)*v_H*G_H(:,1);
s_OH(:,1) = (z_OH^2)*v_OH*G_OH(:,1);
s_C(:,1) = (z_C^2)*v_C*G_C(:,1);
Sigma(:,1) = (F^2)*(s_H(:,1) + s_OH(:,1) + s_C(:,1));
```