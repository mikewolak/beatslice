# Makefile for Beatslice - Audio Beat Detection and Slicing Tool
# 
# Builds with libsndfile and FFTW3 from Homebrew installation

# Project configuration
PROJECT_NAME = beatslice
VERSION = 1.0.0

# Homebrew paths (adjust if needed)
HOMEBREW_PREFIX = /Users/MWOLAK/homebrew
HOMEBREW_INCLUDE = $(HOMEBREW_PREFIX)/include
HOMEBREW_LIB = $(HOMEBREW_PREFIX)/lib

# Compiler and flags
CC = gcc
CFLAGS = -std=c99 -Wall -Wextra -O2 -g
CFLAGS += -I$(HOMEBREW_INCLUDE)
CFLAGS += -DVERSION_MAJOR=1 -DVERSION_MINOR=0 -DVERSION_PATCH=0

# Debug build flags
DEBUG_CFLAGS = -std=c99 -Wall -Wextra -g -O0 -DDEBUG
DEBUG_CFLAGS += -I$(HOMEBREW_INCLUDE)
DEBUG_CFLAGS += -DVERSION_MAJOR=1 -DVERSION_MINOR=0 -DVERSION_PATCH=0

# Linker flags
LDFLAGS = -L$(HOMEBREW_LIB)
LIBS = -lsndfile -lfftw3 -lm

# Installation paths
PREFIX = /usr/local
BINDIR = $(PREFIX)/bin
MANDIR = $(PREFIX)/share/man/man1

# Source files
SRCDIR = src
SOURCES = $(SRCDIR)/main.c \
          $(SRCDIR)/audio_analysis.c \
          $(SRCDIR)/fileio.c \
          $(SRCDIR)/utils.c

# Object files
OBJDIR = obj
OBJECTS = $(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)
DEBUG_OBJECTS = $(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%-debug.o)

# Header files
HEADERS = src/beatslice.h src/audio_analysis.h src/fileio.h src/utils.h

# Target executable
TARGET = $(PROJECT_NAME)
DEBUG_TARGET = $(PROJECT_NAME)-debug

# Test files directory
TESTDIR = test
TEST_AUDIO = $(TESTDIR)/test_audio.wav

# Build rules
.PHONY: all clean debug install uninstall test check-deps help

# Default target
all: check-deps $(TARGET)

# Release build
$(TARGET): $(OBJDIR) $(OBJECTS)
	@echo "Linking $(TARGET)..."
	$(CC) $(OBJECTS) $(LDFLAGS) $(LIBS) -o $(TARGET)
	@echo "Build complete: $(TARGET)"

# Debug build
debug: check-deps $(DEBUG_TARGET)

$(DEBUG_TARGET): $(OBJDIR) $(DEBUG_OBJECTS)
	@echo "Linking $(DEBUG_TARGET)..."
	$(CC) $(DEBUG_OBJECTS) $(LDFLAGS) $(LIBS) -o $(DEBUG_TARGET)
	@echo "Debug build complete: $(DEBUG_TARGET)"

# Object file compilation (release)
$(OBJDIR)/%.o: $(SRCDIR)/%.c $(HEADERS)
	@echo "Compiling $< (release)..."
	$(CC) $(CFLAGS) -c $< -o $@

# Object file compilation (debug)
$(OBJDIR)/%-debug.o: $(SRCDIR)/%.c $(HEADERS)
	@echo "Compiling $< (debug)..."
	$(CC) $(DEBUG_CFLAGS) -c $< -o $@

# Create object directory
$(OBJDIR):
	@mkdir -p $(OBJDIR)

# Check dependencies
check-deps:
	@echo "Checking dependencies..."
	@if [ ! -f "$(HOMEBREW_LIB)/libsndfile.dylib" ] && [ ! -f "$(HOMEBREW_LIB)/libsndfile.so" ]; then \
		echo "Error: libsndfile not found in $(HOMEBREW_LIB)"; \
		echo "Install with: brew install libsndfile"; \
		exit 1; \
	fi
	@if [ ! -f "$(HOMEBREW_LIB)/libfftw3.dylib" ] && [ ! -f "$(HOMEBREW_LIB)/libfftw3.so" ]; then \
		echo "Error: FFTW3 not found in $(HOMEBREW_LIB)"; \
		echo "Install with: brew install fftw"; \
		exit 1; \
	fi
	@if [ ! -f "$(HOMEBREW_INCLUDE)/sndfile.h" ]; then \
		echo "Error: sndfile.h not found in $(HOMEBREW_INCLUDE)"; \
		echo "Install with: brew install libsndfile"; \
		exit 1; \
	fi
	@if [ ! -f "$(HOMEBREW_INCLUDE)/fftw3.h" ]; then \
		echo "Error: fftw3.h not found in $(HOMEBREW_INCLUDE)"; \
		echo "Install with: brew install fftw"; \
		exit 1; \
	fi
	@echo "Dependencies OK"

# Installation
install: $(TARGET)
	@echo "Installing $(TARGET) to $(BINDIR)..."
	install -d $(BINDIR)
	install -m 755 $(TARGET) $(BINDIR)
	@if [ -f "docs/$(PROJECT_NAME).1" ]; then \
		echo "Installing man page..."; \
		install -d $(MANDIR); \
		install -m 644 docs/$(PROJECT_NAME).1 $(MANDIR); \
	fi
	@echo "Installation complete"

# Uninstallation
uninstall:
	@echo "Uninstalling $(TARGET)..."
	rm -f $(BINDIR)/$(TARGET)
	rm -f $(MANDIR)/$(PROJECT_NAME).1
	@echo "Uninstallation complete"

# Create test directory and sample audio
$(TESTDIR):
	@mkdir -p $(TESTDIR)

# Generate test audio (requires SoX)
$(TEST_AUDIO): $(TESTDIR)
	@if command -v sox >/dev/null 2>&1; then \
		echo "Generating test audio file..."; \
		sox -n $(TEST_AUDIO) synth 10 sin 440 vol 0.5; \
	else \
		echo "Warning: SoX not found, cannot generate test audio"; \
		echo "Please place a test audio file at $(TEST_AUDIO)"; \
	fi

# Basic functionality test
test: $(TARGET) $(TEST_AUDIO)
	@echo "Running basic functionality test..."
	@if [ -f "$(TEST_AUDIO)" ]; then \
		echo "Testing with $(TEST_AUDIO)..."; \
		./$(TARGET) --verbose --sensitivity 0.5 $(TEST_AUDIO); \
		echo "Test completed - check output files"; \
	else \
		echo "No test audio file found at $(TEST_AUDIO)"; \
		echo "Run 'make $(TEST_AUDIO)' to generate one, or provide your own"; \
	fi

# Run with test audio and various sensitivity settings
test-sensitivity: $(TARGET) $(TEST_AUDIO)
	@echo "Testing various sensitivity settings..."
	@for sens in 0.1 0.3 0.5 0.7 0.9; do \
		echo "Testing sensitivity $$sens..."; \
		./$(TARGET) --verbose --sensitivity $$sens --output test_sens_$$sens $(TEST_AUDIO); \
	done
	@echo "Sensitivity test complete"

# Memory leak check (requires valgrind)
memcheck: $(DEBUG_TARGET) $(TEST_AUDIO)
	@if command -v valgrind >/dev/null 2>&1; then \
		echo "Running memory leak check..."; \
		valgrind --leak-check=full --show-leak-kinds=all ./$(DEBUG_TARGET) $(TEST_AUDIO); \
	else \
		echo "Valgrind not found, skipping memory check"; \
	fi

# Performance profiling (requires gprof)
profile: CFLAGS += -pg
profile: LDFLAGS += -pg
profile: $(TARGET) $(TEST_AUDIO)
	@echo "Building with profiling support..."
	@$(MAKE) clean
	@$(MAKE) $(TARGET)
	@if [ -f "$(TEST_AUDIO)" ]; then \
		echo "Running profiling test..."; \
		./$(TARGET) $(TEST_AUDIO); \
		gprof $(TARGET) gmon.out > profile_report.txt; \
		echo "Profile report saved to profile_report.txt"; \
	fi

# Static analysis (requires clang-tidy)
analyze:
	@if command -v clang-tidy >/dev/null 2>&1; then \
		echo "Running static analysis..."; \
		clang-tidy $(SOURCES) -- $(CFLAGS); \
	else \
		echo "clang-tidy not found, skipping analysis"; \
	fi

# Code formatting (requires clang-format)
format:
	@if command -v clang-format >/dev/null 2>&1; then \
		echo "Formatting source code..."; \
		clang-format -i $(SOURCES) $(HEADERS); \
		echo "Code formatting complete"; \
	else \
		echo "clang-format not found, skipping formatting"; \
	fi

# Generate documentation (requires doxygen)
docs:
	@if command -v doxygen >/dev/null 2>&1; then \
		echo "Generating documentation..."; \
		doxygen Doxyfile; \
		echo "Documentation generated in docs/html/"; \
	else \
		echo "Doxygen not found, cannot generate documentation"; \
	fi

# Clean build artifacts
clean:
	@echo "Cleaning build artifacts..."
	rm -rf $(OBJDIR)
	rm -f $(TARGET) $(DEBUG_TARGET)
	rm -f *.o core gmon.out profile_report.txt
	rm -f test_*.*  # Remove test output files
	@echo "Clean complete"

# Clean everything including test files
distclean: clean
	@echo "Cleaning all generated files..."
	rm -rf $(TESTDIR)
	rm -rf docs/html docs/latex
	@echo "Distclean complete"

# Show build information
info:
	@echo "Beatslice Build Information"
	@echo "=========================="
	@echo "Project: $(PROJECT_NAME) v$(VERSION)"
	@echo "Compiler: $(CC)"
	@echo "Homebrew prefix: $(HOMEBREW_PREFIX)"
	@echo "Include path: $(HOMEBREW_INCLUDE)"
	@echo "Library path: $(HOMEBREW_LIB)"
	@echo "CFLAGS: $(CFLAGS)"
	@echo "LDFLAGS: $(LDFLAGS)"
	@echo "LIBS: $(LIBS)"
	@echo "Sources: $(SOURCES)"
	@echo "Target: $(TARGET)"

# Help message
help:
	@echo "Beatslice Makefile Help"
	@echo "======================"
	@echo ""
	@echo "Targets:"
	@echo "  all          Build release version (default)"
	@echo "  debug        Build debug version with debug symbols"
	@echo "  install      Install to system (requires sudo)"
	@echo "  uninstall    Remove from system (requires sudo)"
	@echo "  test         Run basic functionality test"
	@echo "  test-sensitivity Test various sensitivity settings"
	@echo "  memcheck     Run memory leak check (requires valgrind)"
	@echo "  profile      Build with profiling and generate report"
	@echo "  analyze      Run static code analysis (requires clang-tidy)"
	@echo "  format       Format source code (requires clang-format)"
	@echo "  docs         Generate documentation (requires doxygen)"
	@echo "  clean        Remove build artifacts"
	@echo "  distclean    Remove all generated files"
	@echo "  info         Show build configuration"
	@echo "  help         Show this help message"
	@echo ""
	@echo "Prerequisites:"
	@echo "  brew install libsndfile fftw"
	@echo ""
	@echo "Example usage:"
	@echo "  make                    # Build release version"
	@echo "  make debug              # Build debug version"
	@echo "  make test               # Run basic test"
	@echo "  sudo make install       # Install system-wide"

# Packaging (create distributable archive)
PACKAGE_NAME = $(PROJECT_NAME)-$(VERSION)
package: clean
	@echo "Creating package $(PACKAGE_NAME)..."
	@mkdir -p $(PACKAGE_NAME)
	@cp -r src/ $(PACKAGE_NAME)/
	@cp Makefile $(PACKAGE_NAME)/
	@cp README.md $(PACKAGE_NAME)/ 2>/dev/null || true
	@cp LICENSE $(PACKAGE_NAME)/ 2>/dev/null || true
	@tar czf $(PACKAGE_NAME).tar.gz $(PACKAGE_NAME)
	@rm -rf $(PACKAGE_NAME)
	@echo "Package created: $(PACKAGE_NAME).tar.gz"
