% Function calculates Heat flow from Human Body to TEG Surface

 function Q = ContactRes(Tskin, Thg)

 global hxlC hxwC
 
if Tskin < Thg
    fprintf('Th Guess too big\n')
    return
else
end



% hxlC = 0.04;
% hxwC = 0.04;


Area = hxlC*hxwC;

hcontact = 16.667;

Rth = 1/(Area*hcontact);

Q = (Tskin - Thg)/Rth;

