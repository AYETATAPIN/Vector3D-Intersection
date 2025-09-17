#include <iostream>
#include <cmath>
#include <string>
#include <stdexcept>
#include <sstream>
#include <gtest/gtest.h>

using namespace std;

constexpr double INACCURACY = 1e-10;

class Vector3D {
private:
    double X;
    double Y;
    double Z;

public:
    Vector3D(const double x, const double y, const double z) {
        this->X = x;
        this->Y = y;
        this->Z = z;
    }

    void setX(const double X) noexcept {
        this->X = X;
    }

    void setY(const double Y) noexcept {
        this->Y = Y;
    }

    void setZ(const double Z) noexcept {
        this->Z = Z;
    }

    double getX() const noexcept {
        return this->X;
    }

    double getY() const noexcept {
        return this->Y;
    }

    double getZ() const noexcept {
        return this->Z;
    }

    Vector3D operator+(const Vector3D &vector_3d) const noexcept {
        return Vector3D(
            this->X + vector_3d.getX(),
            this->Y + vector_3d.getY(),
            this->Z + vector_3d.getZ()
        );
    }

    Vector3D operator-(const Vector3D &vector_3d) const noexcept {
        return Vector3D(
            this->X - vector_3d.getX(),
            this->Y - vector_3d.getY(),
            this->Z - vector_3d.getZ()
        );
    }

    Vector3D operator*(const double scalar) const noexcept {
        return Vector3D(
            this->X * scalar,
            this->Y * scalar,
            this->Z * scalar
        );
    }

    bool operator==(const Vector3D &vector_3d) const noexcept {
        return fabs(this->X - vector_3d.X) < INACCURACY
               && fabs(this->Y - vector_3d.Y) < INACCURACY
               && fabs(this->Z - vector_3d.Z) < INACCURACY;
    }

    bool operator!=(const Vector3D &vector_3d) const noexcept {
        return !(*this == vector_3d);
    }

    explicit operator string() const noexcept {
        return "{" + to_string(this->X) + ", " + to_string(this->Y) + ", " + to_string(this->Z) + "}";
    }

    double length() const noexcept {
        return sqrt(X * X + Y * Y + Z * Z);
    }

    bool is_uniform() const noexcept {
        return fabs(this->X - this->Y) < INACCURACY && fabs(this->Y - this->Z) < INACCURACY;
    }
};

Vector3D vector_product(const Vector3D &first, const Vector3D &second) noexcept {
    return Vector3D(
        first.getY() * second.getZ() - first.getZ() * second.getY(),
        first.getZ() * second.getX() - first.getX() * second.getZ(),
        first.getX() * second.getY() - first.getY() * second.getX()
    );
}

double scalar_product(const Vector3D &first, const Vector3D &second) noexcept {
    return first.getX() * second.getX() + first.getY() * second.getY() + first.getZ() * second.getZ();
}

class Segment3D {
private:
    Vector3D start;
    Vector3D end;

public:
    Segment3D(const Vector3D &start, const Vector3D &end): start(start), end(end) {
    }

    void setStart(const Vector3D &start) noexcept {
        this->start = start;
    }

    void setEnd(const Vector3D &end) noexcept {
        this->end = end;
    }

    Vector3D getStart() const noexcept {
        return this->start;
    }

    Vector3D getEnd() const noexcept {
        return this->end;
    }

    bool operator==(const Segment3D &segment_3d) const noexcept {
        return this->start == segment_3d.getStart() && this->end == segment_3d.getEnd();
    }

    Vector3D directional_vector() const noexcept {
        return end - start;
    }

    double length() const noexcept {
        return this->directional_vector().length();
    }
};

Vector3D Intersect(const Segment3D &first, const Segment3D &second) {
    if (first == second) {
        throw invalid_argument(string("Segments are the same and the intersection is an interval with start ")
                               + string(first.getStart())
                               + string(" and end ")
                               + string(second.getEnd()));
    }

    const Vector3D start_diff = second.getStart() - first.getStart();
    const Vector3D end_diff = second.getEnd() - first.getEnd();

    // Collinear check
    if (start_diff.is_uniform() && end_diff.is_uniform()) {
        if (first.getEnd().getX() < second.getStart().getX() || first.getStart().getX() > second.getEnd().getX()) {
            throw invalid_argument(string("Segments are on the same line and do not intersect"));
        }

        const double interval_x1 = max(first.getStart().getX(), second.getStart().getX());
        const double interval_y1 = max(first.getStart().getY(), second.getStart().getY());
        const double interval_z1 = max(first.getStart().getZ(), second.getStart().getZ());

        const double interval_x2 = min(first.getEnd().getX(), second.getEnd().getX());
        const double interval_y2 = min(first.getEnd().getY(), second.getEnd().getY());
        const double interval_z2 = min(first.getEnd().getZ(), second.getEnd().getZ());

        const Vector3D startInterval(interval_x1, interval_y1, interval_z1);
        const Vector3D endInterval(interval_x2, interval_y2, interval_z2);

        throw invalid_argument(string("Segments intersection is an interval with start ")
                               + string(startInterval)
                               + string(" and end ")
                               + string(endInterval));
    }

    if (first.getStart() == second.getStart() || first.getStart() == second.getEnd()) {
        return first.getStart();
    }
    if (first.getEnd() == second.getEnd() || first.getEnd() == second.getStart()) {
        return first.getEnd();
    }

    const Vector3D first_direction = first.directional_vector();
    const Vector3D second_direction = second.directional_vector();

    const Vector3D directions_vector_product = vector_product(first_direction, second_direction);

    if (directions_vector_product.length() < INACCURACY) {
        throw invalid_argument(string("Segments are parallel and do not intersect"));
    }

    const Vector3D cross1 = vector_product(start_diff, second_direction);
    const Vector3D cross2 = vector_product(start_diff, first_direction);

    const double t1 = scalar_product(cross1, directions_vector_product) / scalar_product(
                          directions_vector_product, directions_vector_product);
    const double t2 = scalar_product(cross2, directions_vector_product) / scalar_product(
                          directions_vector_product, directions_vector_product);

    if (t1 >= 0 && t1 <= 1 && t2 >= 0 && t2 <= 1) {
        return first.getStart() + first_direction * t1;
    }

    throw invalid_argument(string("Segments do not intersect"));
}


TEST(VectorOperations, VectorProduct) {
    const Vector3D first(1, 0, 0);
    const Vector3D second(0, 1, 0);
    EXPECT_EQ(vector_product(first, second), Vector3D(0, 0, 1));
}

TEST(VectorOperations, ScalarProduct) {
    const Vector3D first(1, 2, 3);
    const Vector3D second(4, 5, 6);
    EXPECT_NEAR(scalar_product(first, second), 32, INACCURACY);
}

TEST(SegmentIntersection, BasicIntersection) {
    const Vector3D start1(0, 0, 0);
    const Vector3D end1(1, 1, 0);
    const Vector3D start2(0, 1, 0);
    const Vector3D end2(1, 0, 0);

    const Segment3D segment1(start1, end1);
    const Segment3D segment2(start2, end2);

    const Vector3D intersection = Intersect(segment1, segment2);
    EXPECT_EQ(intersection, Vector3D(0.5, 0.5, 0));
}

TEST(SegmentIntersection, CollinearSegments) {
    const Vector3D start1(0, 0, 0);
    const Vector3D end1(1, 0, 0);
    const Vector3D start2(0, 1, 0);
    const Vector3D end2(1, 1, 0);
    const Segment3D segment1(start1, end1);
    const Segment3D segment2(start2, end2);

    EXPECT_THROW({
                 Intersect(segment1, segment2);
                 }, invalid_argument);
}

TEST(SegmentIntersection, SameSegments) {
    const Vector3D start1(0, 0, 0);
    const Vector3D end1(1, 0, 0);
    const Vector3D start2(0, 0, 0);
    const Vector3D end2(1, 0, 0);
    const Segment3D segment1(start1, end1);
    const Segment3D segment2(start2, end2);

    EXPECT_THROW({
                 Intersect(segment1, segment2);
                 }, invalid_argument);
}

TEST(SegmentIntersection, SameLineNoIntersection) {
    const Vector3D start1(0, 0, 0);
    const Vector3D end1(1, 1, 1);
    const Vector3D start2(2, 2, 2);
    const Vector3D end2(3, 3, 3);
    const Segment3D segment1(start1, end1);
    const Segment3D segment2(start2, end2);

    EXPECT_THROW({
                 Intersect(segment1, segment2);
                 }, invalid_argument);
}

TEST(SegmentIntersection, SameLineWithIntersection) {
    const Vector3D start1(0, 0, 0);
    const Vector3D end1(1, 1, 1);
    const Vector3D start2(0.5, 0.5, 0.5);
    const Vector3D end2(3, 3, 3);
    const Segment3D segment1(start1, end1);
    const Segment3D segment2(start2, end2);

    EXPECT_THROW({
                 Intersect(segment1, segment2);
                 }, invalid_argument);
}

TEST(SegmentIntersection, NoIntersection) {
    const Vector3D start1(0, 0, 0);
    const Vector3D end1(1, 1, 1);
    const Vector3D start2(227, 228, 737);
    const Vector3D end2(1337, 911, 420);
    const Segment3D segment1(start1, end1);
    const Segment3D segment2(start2, end2);

    EXPECT_THROW({
                 Intersect(segment1, segment2);
                 }, invalid_argument);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
