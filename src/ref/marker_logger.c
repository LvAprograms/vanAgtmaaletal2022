#define XINTERVAL 50e3

#ifndef STM_MAX_MARKER_LOG_CAPACITY
#define STM_MAX_MARKER_LOG_CAPACITY 50
#endif

static int const marker_log_capacity = STM_MAX_MARKER_LOG_CAPACITY; // defined twice? maybe remove later
static struct MarkerLogConfig marker_log_configs[STM_MAX_MARKER_LOG_CAPACITY]; // limit how many configs are initialised

/*
---- RECTANGLE RELATED STRUCTURES/FUNCTIONS ----
*/
struct Rectangle {
  double xmin, xmax, ymin, ymax;
};

int marker_in_rectangle(struct Rectangle *const rect, int const index) {
	// checks if a marker is within the predefined rectangle 
  int ret = 0;
  double const x = markx[ii], y = marky[ii];
  if (x > rect->xmin && x < rect->xmax && y > rect->ymin && y < rect->ymax) {
    ret = 1;
  }
  return ret;
}

enum Found { Found, NotFound };

/* 
---- LINE READING RELATED FUNCTIONS ----
*/
const char *read_line(FILE *file) {
  // reads a line from the input file 
  
  // allocate memory
  int maximumLineLength = 128;
  char *lineBuffer = (char *)malloc(sizeof(char) * maximumLineLength);
  if (lineBuffer == NULL) {
    printf("Error allocating memory for line buffer.");
    exit(1);
  }

  char ch = getc(file); // read the first character of the line
  int count = 0;


  while ((ch != '\n') && (ch != EOF)) {
    if (count == maximumLineLength) {
      maximumLineLength += 128; // why would you increase the buffer size?
      lineBuffer = realloc(lineBuffer, maximumLineLength);
      if (lineBuffer == NULL) {
        printf("Error reallocating space for line buffer.");
        exit(1);
      }
    }
    lineBuffer[count] = ch;
    count++;
    ch = getc(file);
  }

  lineBuffer[count] = '\0'; //make sure the string has the null-pointer?
  return lineBuffer;
}

int is_white_space(char const c) { return c == ' ' || c == '\t' || c == '\n'; };

int is_comment(char const *const line) {
  int const len = strlen(line);
  for (int ii = 0; ii < len; ii++) {
    if (is_white_space(line[ii])) {
      continue;
    }
    if (line[ii] == '#') {
      return 1;
    }
  }
  return 0;
}
/* 
---- END LINE READING RELATED FUNCTIONS ----
---- MARKER LOGGING RELATED FUNCTIONS ----
*/
struct MarkerLogConfig {
  // configuration structure: contains rock type and a rectangle
  int rock_id;
  struct Rectangle rect;
};



static void check_marker_log_capacity(int const curr) {
  if (curr > marker_log_capacity) {
    printf("Maximum marker log capacity (%d) exceeded with: %d\n",
           marker_log_capacity, curr);
    exit(1);
  }
}

// file structure of a marker config file to be loaded
/*
 * # rock_id, xmin, xmax, ymin, ymax
 */
int load_marker_log_config(FILE *file) {
  int n_read = 0; // the number of markers read

  while (1) {
    if (feof(file)) {
      break;
    }

    char const *line = read_line(file); // read a line from the config file
    if (NULL == line) {
      return -1;
    }

    if (is_comment(line)) {
      free(line);
      continue;
    }

	// put the line buffer into the MarkerLogConfig:
    struct MarkerLogConfig current = {};
    struct Rectangle *const rect = &current.rect;
    int const ret = sscanf(line, "%d %lf %lf %lf %lf", &current.rock_id,
                           &rect->xmin, &rect->xmax, &rect->ymin, &rect->ymax);

    if (ret <= 0) {
      printf("error processing line '%s'\n", line);
      free(line);
      exit(1);
    }
    check_marker_log_capacity(n_read);

    marker_log_configs[n_read] = current;
    n_read += 1;

    free(line);
  }

  return n_read;
}

int search_for_marker(struct MarkerLogConfig *const cfg) {
  for (int ii = 0; ii < marknum; ii++) {
    if ((cfg->rock_id == markt[ii]) && marker_in_rectangle(&cfg->rect, ii)) {
      return ii;
    }
  }

  return -1;
}

void print_found_markers(int const n_read) {
  printf("Searching for markers...\n");
  int n_found = 0;
  for (int ii = 0; ii < n_read; ii++) {
    int const index = search_for_marker(&marker_log_configs[ii]);
    if (index >= 0) {
      n_found += 1;
      printf("%d\n", index);
    }
  }
  printf("Search complete, found: %d\n", n_found);
}

int search_markers(int argc, char **argv) {
  FILE *config_file = fopen();
  // check for NULL

  int const n_read = load_marker_log_config(config_file);
  print_found_markers(n_read);

  fclose(config_file);

  return 0;
}

static int marker_indices[STM_MAX_MARKER_LOG_CAPACITY];

int load_marker_indices_from_file(FILE *file) {
  int n_read = 0;
  while (1) {
    check_marker_log_capacity(n_read);

    int const ret = fscanf(file, "%d\n", &marker_indices[n_read]);
    if (ret <= 0) {
      printf("failed to read marker index from file.\n");
      exit(1);
    }
    n_read += 1;
  }
  return n_read;
}

void print_marker_data_file(FILE *file, int const n_read) {
  for (int ii = 0; ii < n_read; ii++) {
    int const index = marker_indices[ii];
    fprintf(file, "");
  }
}

int print_marker_data(int argc, char **argv) {
  FILE *marker_indices = fopen(argv[0]); 
  if (marker_indices == NULL) {
	printf("error opening marker indices file %s\n", argv[0]);
	exit(1);
  }

  int const n_read = load_marker_indices_from_file(marker_indices);
  fclose(marker_indices);

  FILE *out = fopen("", "a");
  print_marker_data_file(stdout, n_read);

  fclose(out);
  return 0;
}

struct MarkerLogger {
  enum Found is_found;
  int indices[5];
  int rock_id;
  struct Rectangle rect;
};

typedef int (*Run)(int, char **); // function pointer

int main(int argc, char** argv) {
	char* command_names[][] = {
		"search_markers",
		"print_marker_data",
	};
	Run commands[10] = {
		search_markers,
		print_marker_data,
	};
	// ./main search _markers
	// 0.  1.
	for (int ii == 0; ii < 2; ii++) {
	if (strcmp(argv[1], command_names[ii]) == 0) {
		return commands[ii](argc, argv);
	}
	}
}