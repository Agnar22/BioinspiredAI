#ifndef DIR_H
#define DIR_H

enum class Dir:int {
    s=0, u, l, d, r
};

static Dir find_direction(int from, int to, int width, int height) {
    if (!(0 <= from && from < width*height && 0 <= to && to < width*height))
        throw std::invalid_argument("The given points were outside of the valid range.");
    if (to-from == 0)
        return Dir::s;
    if (std::abs(to-from) == width)
        return to>from? Dir::d : Dir::u;
    if (to-from == 1 && from%width!=width-1)
        return Dir::r;
    if (to-from == -1 && from%width!=0)
        return Dir::l;
    throw std::invalid_argument("There were no orthogonal direction between the two points.");
};

#endif