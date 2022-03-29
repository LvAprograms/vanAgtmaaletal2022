#define XINTERVAL 50e3
enum IsSet {
    Set,
    NotSet
};

struct Rectangle {
    enum IsSet is_set;
    double xmin, xmax, ymin, ymax;
};

enum Found {
    Found,
    NotFound
};

struct MarkerLogger {
    enum Found is_found;
    int indices[5];
    int rock_id;
    struct Rectangle rect;
};

int search_for_marker(struct MarkerLogger*const logger, int index) {
    if (logger->is_found == Found) {
        return logger->indices[index];
    }

    for (int ii = 0; ii < marknum; ii++) {
        if (logger->rock_id != markt[ii]) {
            continue;
        }

        struct Rectangle const*const rect = &logger->rect;
        if (rect->is_set == NotSet) {
            continue;
        }

        double const x = markx[ii], y = marky[ii];
        if (x > (rect->xmin + index * XINTERVAL) && x < rect->xmax &&
            y > rect->ymin && y < rect->ymax) {
            if (index == 4) logger->is_found = Found;
            logger->indices[index] = ii;
            return ii;
        }
    }

    return -1;
}

/*
Problem: the current implementation only logs ONE marker. I want to log multiple, with a step in X direction
*/
