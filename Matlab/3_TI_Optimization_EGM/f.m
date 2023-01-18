function res = f(capital,wealth)

global beta gamma alpha delta;

if capital > 0
    res = capital.^alpha + (1-delta)*capital - wealth;
else
   res = 100000;
end