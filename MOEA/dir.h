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

static int find_pos(int from_pos, Dir dir, int width, int height) {
    if (!(0 <= from_pos && from_pos < width*height))
        throw std::invalid_argument("The given point is outside of the valid range.");
    int x = from_pos%width;
    int y = from_pos/width;
    if (x==0 && dir==Dir::l || x==width-1 && dir==Dir::r ||y==0 && dir==Dir::u || y==height-1 && dir==Dir::d)
        throw std::invalid_argument("The direction was not possible from the given point.");
    if (dir==Dir::s)
        return from_pos;
    if (dir==Dir::u)
        return from_pos - width;
    if (dir==Dir::d)
        return from_pos + width;
    if (dir==Dir::l)
        return from_pos - 1;
    if (dir==Dir::r)
        return from_pos + 1;
    throw std::invalid_argument("No direction was found.");
}

static Dir reverse_dir(Dir dir) {
    if (dir==Dir::s)
        return Dir::s;
    if (dir==Dir::u)
        return Dir::d;
    if (dir==Dir::d)
        return Dir::u;
    if (dir==Dir::l)
        return Dir::r;
    if (dir==Dir::r)
        return Dir::l;
    throw std::invalid_argument("No direction was found.");
}

#endif