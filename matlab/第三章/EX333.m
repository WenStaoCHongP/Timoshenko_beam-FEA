mode=1;
for loopi=1:elment_number
    x=[nc(en(loopi,1),1) nc(en(loopi,2),1)];
    y=[nc(en(loopi,1),2) nc(en(loopi,2),2)];
    line(x,y)
    mx(1:2)=0;
    my(1:2)=0;
    for loopj=1:2
        if ed(en(loopi,loopj),1)~=0
            mx(loopj)=x(loopj)+v(ed(en(loopi,loopj),1),mode)/8;
        else
            mx(loopj)=x(loopj);
        end
         if ed(en(loopi,loopj),2)~=0
            my(loopj)=y(loopj)+v(ed(en(loopi,loopj),2),mode)/8;
        else
            my(loopj)=y(loopj);
         end
    end
    line(mx,my,'LineStyle',':')
end
