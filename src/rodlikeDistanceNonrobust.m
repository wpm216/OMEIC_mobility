function distance = rodlikeDistanceNonrobust(p, q)
% Taken directly from 
%   https://www.geometrictools.com/Documentation/DistanceLine3Line3.pdf.
% Find the distance between two rods (3D).
% p and q are 2x3 arrays with the beginning and end points of each segment.
% This algorithm is "not robust", meaning it may fail in some cases 
%  (see Examples 1-2, pp.7-8). No worrying cases 
%  (of closely parallel segments) were detected in our tests, so this 
%  algorithm will suffice.

p0 = p(1:3);
p1 = p(4:6);
q0 = q(1:3);
q1 = q(4:6);
a = (p1-p0) * (p1-p0).';
b = (p1-p0) * (q1-q0).';
c = (q1-q0) * (q1-q0).';
d = (p1-p0) * (p0-q0).';
e = (q1-q0) * (p0-q0).';
% f = (p0-q0) * (p0-q0).';

den = a*c - b^2;

if den > 10^-8 % line segments are not approximately parallel
    bte = b*e; ctd = c*d;
    if bte <= ctd % s <= 0
        if e <= 0 % t <= 0 (region 6)
            if -d >= a
                s = 1;
            elseif -d > 0
                s = -d/a;
            else
                s = 0;
            end
            t = 0;
        elseif e < c % 0 < t < 1 (region 5)
            s = 0;
            t = e/c;
        else % t >=1 (region 4)
            if b-d >= a
                s = 1;
            elseif b-d > 0
                s = (b-d)/a;
            else
                s = 0;
            end
            t = 1;
        end
    else % s < 0
        s = bte - ctd;
        if s >= den % s >= 1
            if b+e <= 0 % t <= 0 (region 8)
                if -d <= 0
                    s = 0;
                elseif -d < a
                    s = -d/a;
                else
                    s = 1;
                end
                t = 0;
            elseif b+e < c % 0 < t < 1 (region 1)
                s = 1;
                t = (b+e)/c;
            else % t >= 1 (region 2)
                if b-d <= 0
                    s = 0;
                elseif b-d < a
                    s = (b-d)/a;
                else
                    s = 1;
                end
                t = 1;
            end
        else % 0 < s < 1
            ate = a*e; btd = b*d;
            if ate <= btd % t <= 0 (region 7)
                if -d <= 0
                    s = 0;
                elseif -d >= a
                    s = 1;
                else
                    s = -d/a;
                end
                t = 0;
            else % t > 0
                t = ate-btd;
                if t >= den % t >=1 (region 3)
                    if b-d <= 0
                        s = 0;
                    elseif b-d >= a
                        s = 1;
                    else
                        s = (b-d)/a;
                    end
                    t = 1;
                else % 0 < t < 1 (region 0)
                    s = s/den;
                    t = t/den;
                end
            end
        end
    end

else % parallel segments,
     % (https://www.geometrictools.com/Documentation/DistanceLine3Line3.pdf
     % pp.6-7)
    if e <= 0
        t = 0;
        if -d <= 0
            s = 0;
        elseif -d >= a
            s = 1;
        else
            s = -d/a;
        end
    elseif e >= c
        t = 1;
        if b-d <= 0
            s = 0;
        elseif b-d >= a
            s = 1;
        else
            s = (b-d)/a;
        end
    else 
        s = 0;
        t = e/c;
    end 
end

p2 = (1-s)*p0 + s*p1;
q2 = (1-t)*q0 + t*q1;
distance = norm(p2-q2);




    