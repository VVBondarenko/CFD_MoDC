#ifndef STREAMLINE_H
#define STREAMLINE_H

class StreamLinePoint
{
public:
    StreamLinePoint(double x, double y);

    double x,y;
    void move(double dx, double dy);
};

#endif // STREAMLINE_H
