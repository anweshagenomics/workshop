# Unix & Linux Basics — Day 1

> **Five-day Hands-on Workshop Program on Bioinformatics for One Health and Pathogen-Specific Research**  
> **Instructor:** Anwesha  
> **Duration:** 3 hours (180 minutes)  
> **Level:** Beginner  
> **Organized by:** CSIR–Institute of Genomics & Integrative Biology (IGIB)

A beginner's guide to the command line, file systems, and essential Unix/Linux concepts — all commands run directly in your terminal.

---

## Table of Contents

1. [Origins — Unix and Linux](#1-origins--unix-and-linux)
2. [Key Terminologies](#2-key-terminologies)
3. [Why the Command Line?](#3-why-the-command-line)
4. [Navigation Commands](#4-navigation-commands)
5. [File System & Paths](#5-file-system--paths)
6. [Tab Completion](#6-tab-completion)
7. [Creating Directories](#7-creating-directories)
8. [Creating Files](#8-creating-files)
9. [Copying and Moving Files](#9-copying-and-moving-files)
10. [Removing Files and Directories](#10-removing-files-and-directories)
11. [Wildcards](#11-wildcards)
12. [Viewing File Contents](#12-viewing-file-contents)
13. [Counting with wc](#13-counting-with-wc)
14. [Output Redirection](#14-output-redirection)
15. [Pipes](#15-pipes)
16. [Working with CSV Data](#16-working-with-csv-data)
17. [Command Chaining](#17-command-chaining)
18. [Best Practices](#18-best-practices)
19. [Practice Exercises](#19-practice-exercises)
20. [Summary](#20-summary)
21. [Commands Reference Sheet](#commands-reference)
---

## 1. Origins — Unix and Linux

| Era | Event |
|-----|-------|
| 1970s | Unix created at Bell Labs by Ken Thompson & Dennis Ritchie |
| 1991 | Linux kernel released by Linus Torvalds as open-source |
| Today | Powers smartphones, supercomputers, and cloud infrastructure |

---

## 2. Key Terminologies

| Term | What it is |
|------|------------|
| **Unix** | Family of multitasking, multi-user operating systems (Linux, macOS) |
| **Terminal** | The application where you type commands (e.g. GNOME Terminal, iTerm2) |
| **Shell** | The program that interprets your commands (sits between you and the OS) |
| **Bash** | *Bourne Again SHell* — the most common shell, default on most Linux systems |

---

## 3. Why the Command Line?

GUIs are fine for simple tasks, but the CLI is essential when you need to:

- Process thousands of files in seconds
- Automate repetitive tasks with scripts
- Work on remote servers and HPC clusters
- Build reproducible bioinformatics pipelines

> **Case sensitivity matters:** `File` and `file` are two different things in Unix.

---

## 4. Navigation Commands

### Syntax

```
command -options argument
```

**Example:** `ls -l Documents` — list contents of `Documents` in long format

### Commands

| Command | Description |
|---------|-------------|
| `pwd` | Print current working directory |
| `ls` | List files and directories |
| `ls -l` | Detailed (long) listing |
| `ls -a` | Show all files, including hidden ones |
| `ls -la` | Long listing + hidden files |
| `cd <dir>` | Change to a directory |
| `cd ..` | Go up one level |
| `cd ~` | Go to your home directory |
| `cd /` | Go to the root directory |

### Try it

```bash
# Where am I?
pwd

# What's here?
ls -l

# Show hidden files too
ls -la

# Go home and look around
cd ~ && pwd && ls

# Get help for any command
ls --help
man ls
```

---

## 5. File System & Paths

The filesystem organizes data as a tree, rooted at `/`.

```
/
├── home/
│   └── ubuntu/
│       ├── Documents/
│       └── Downloads/
├── etc/
└── usr/
```

### Absolute vs Relative Paths

```bash
# Absolute path — starts from root
/home/ubuntu/Documents/data.txt

# Relative path — relative to where you are now
Documents/data.txt
```

### Special Shortcuts

| Symbol | Meaning |
|--------|---------|
| `~` | Your home directory (e.g. `/home/ubuntu`) |
| `..` | Parent directory (one level up) |
| `.` | Current directory |

---

## 6. Tab Completion

Press **Tab** while typing to auto-complete filenames and commands.

- **One match** → completes automatically
- **Multiple matches** → press Tab twice to see all options

```bash
# Type this and press Tab:
ls Doc<TAB>
# Becomes:
ls Documents/
```

> Make this a habit — it prevents typos and saves enormous time.

---

## 7. Creating Directories

| Command | Description |
|---------|-------------|
| `mkdir <dir>` | Create a new directory |
| `mkdir -p <path>` | Create directory and any missing parents |

```bash
# Create our working directory
mkdir -p ~/bioinfo_day1
cd ~/bioinfo_day1 && pwd

# Create a subdirectory
mkdir ~/bioinfo_day1/thesis_notes
ls -l ~/bioinfo_day1/

# Create nested directories in one command
mkdir -p ~/bioinfo_day1/projects/experiment1
ls -l ~/bioinfo_day1/
```

---

## 8. Creating Files

| Command | Description |
|---------|-------------|
| `touch <file>` | Create an empty file (or update timestamp) |
| `echo "text" > file` | Write text to a file |
| `echo "text" >> file` | Append text to a file |
| `cat > file << 'EOF' ... EOF` | Write a multi-line file (heredoc) |

```bash
# Create multiple empty files at once
cd ~/bioinfo_day1/thesis_notes
touch book1.txt book2.txt notes.txt references.txt
ls -l

# Write content to a file
echo "Chapter 1: Introduction" > book1.txt
echo "Chapter 2: Methods" >> book1.txt
cat book1.txt

# Write a multi-line file using heredoc
cat > sample_data.txt << 'EOF'
Line 1: Introduction to Unix
Line 2: Understanding File Systems
Line 3: Command Line Basics
Line 4: File Management
Line 5: Text Processing
Line 6: Advanced Scripting
Line 7: System Administration
Line 8: Networking Concepts
Line 9: Security Best Practices
Line 10: Advanced Topics
EOF

cat sample_data.txt
```

---

## 9. Copying and Moving Files

| Command | Description |
|---------|-------------|
| `cp <src> <dest>` | Copy a file |
| `cp -r <src> <dest>` | Copy a directory recursively |
| `mv <old> <new>` | Move or rename a file/directory |

```bash
# Copy a file
cp ~/bioinfo_day1/thesis_notes/book1.txt \
   ~/bioinfo_day1/thesis_notes/book1_backup.txt
ls ~/bioinfo_day1/thesis_notes/

# Copy an entire directory
cp -r ~/bioinfo_day1/thesis_notes \
      ~/bioinfo_day1/thesis_notes_backup
ls -l ~/bioinfo_day1/

# Rename a file with mv
mv ~/bioinfo_day1/thesis_notes/book1_backup.txt \
   ~/bioinfo_day1/thesis_notes/archive_book1.txt
ls ~/bioinfo_day1/thesis_notes/
```

---

## 10. Removing Files and Directories

| Command | Description |
|---------|-------------|
| `rm <file>` | Remove a file |
| `rm -i <file>` | Remove with confirmation prompt |
| `rm -r <dir>` | Remove a directory and all its contents |
| `rm -r -i <dir>` | Recursive removal with confirmation |
| `rmdir <dir>` | Remove an **empty** directory only |

> ⚠️ **Warning:** The Unix shell has no Recycle Bin. Deleted files are gone immediately.  
> Use `rm -i` or `rm -r -i` when in doubt.

```bash
# Remove a single file
rm ~/bioinfo_day1/thesis_notes/archive_book1.txt
ls ~/bioinfo_day1/thesis_notes/

# Remove a directory
rm -r ~/bioinfo_day1/thesis_notes_backup
ls ~/bioinfo_day1/
```

---

## 11. Wildcards

Wildcards let you match multiple files at once.

### `*` — matches zero or more characters

| Pattern | Matches |
|---------|---------|
| `*.pdb` | All files ending in `.pdb` |
| `p*.pdb` | Files starting with `p`, ending in `.pdb` |
| `data*` | All files starting with `data` |

### `?` — matches exactly one character

| Pattern | Matches |
|---------|---------|
| `?ethane.pdb` | `methane.pdb` (one char before `ethane`) |
| `???ane.pdb` | `cubane.pdb`, `ethane.pdb`, `octane.pdb` |
| `sample?.txt` | `sample1.txt`, `sample2.txt`, etc. |

```bash
# Create sample files
cd ~/bioinfo_day1
touch methane.pdb ethane.pdb propane.pdb cubane.pdb octane.pdb

# Match all .pdb files
ls *.pdb

# Match only files starting with 'p'
ls p*.pdb

# Create and test ? wildcard
touch sample1.txt sample2.txt sample3.txt data_a.txt data_bb.txt
echo "Files matching sample?.txt:"
ls sample?.txt
```

---

## 12. Viewing File Contents

| Command | Description |
|---------|-------------|
| `cat <file>` | Print entire file to terminal |
| `head <file>` | Show first 10 lines |
| `head -n N <file>` | Show first N lines |
| `tail <file>` | Show last 10 lines |
| `tail -n N <file>` | Show last N lines |
| `less <file>` | Paginated viewer (use arrow keys, `q` to quit) |
| `more <file>` | Simpler paginated viewer |
| `grep <pattern> <file>` | Print lines matching a pattern |

```bash
# Print the whole file
cat ~/bioinfo_day1/thesis_notes/sample_data.txt

# First 3 lines
echo "First 3 lines:"
head -n 3 ~/bioinfo_day1/thesis_notes/sample_data.txt

# Last 3 lines
echo "Last 3 lines:"
tail -n 3 ~/bioinfo_day1/thesis_notes/sample_data.txt

# Search for a pattern
echo "Lines containing 'File':"
grep "File" ~/bioinfo_day1/thesis_notes/sample_data.txt

# Interactive viewer (best for large files)
less ~/bioinfo_day1/thesis_notes/sample_data.txt
```

---

## 13. Counting with `wc`

`wc` — word count. Output format: `lines  words  characters  filename`

| Option | Description |
|--------|-------------|
| `wc <file>` | Lines, words, and characters |
| `wc -l` | Lines only |
| `wc -w` | Words only |
| `wc -c` | Characters (bytes) only |

```bash
# Full count
echo "Complete word count:"
wc ~/bioinfo_day1/thesis_notes/sample_data.txt

# Lines only
wc -l ~/bioinfo_day1/thesis_notes/sample_data.txt

# Words only
wc -w ~/bioinfo_day1/thesis_notes/sample_data.txt

# Characters only
wc -c ~/bioinfo_day1/thesis_notes/sample_data.txt

# Count across multiple files with wildcard
echo "Line count for all .pdb files:"
wc -l ~/bioinfo_day1/*.pdb
```

---

## 14. Output Redirection

| Operator | Description |
|----------|-------------|
| `>` | Write output to a file (**overwrites** existing content) |
| `>>` | **Append** output to a file |
| `<` | Read input from a file |

```bash
# Save ls output to a file
ls ~/bioinfo_day1 > ~/bioinfo_day1/file_list.txt
cat ~/bioinfo_day1/file_list.txt

# Save wc output to a file
wc -l ~/bioinfo_day1/*.pdb > ~/bioinfo_day1/line_counts.txt
cat ~/bioinfo_day1/line_counts.txt

# Append txt file counts (does NOT overwrite)
wc -l ~/bioinfo_day1/*.txt >> ~/bioinfo_day1/line_counts.txt
cat ~/bioinfo_day1/line_counts.txt
```

---

## 15. Pipes

The pipe operator `|` passes the output of one command as input to the next.

```bash
command1 | command2 | command3
```

```bash
# Count how many .pdb files exist
ls ~/bioinfo_day1/*.pdb | wc -l

# Sort files by line count
wc -l ~/bioinfo_day1/*.txt | sort

# Most common pattern: grep + wc
grep "Line" ~/bioinfo_day1/thesis_notes/sample_data.txt | wc -l
```

> The Unix philosophy: small tools that do one thing well, chained together into powerful pipelines.

---

## 16. Working with CSV Data

A typical bioinformatics workflow: create → combine → filter → count → save.

```bash
# Create sample variant CSV files
cat > ~/bioinfo_day1/variants1.csv << 'EOF'
Sample,Variant,Country
Sample1,Alpha,USA
Sample2,Beta,UK
Sample3,Alpha,India
EOF

cat > ~/bioinfo_day1/variants2.csv << 'EOF'
Sample4,Gamma,Brazil
Sample5,Alpha,Canada
Sample6,Delta,India
EOF

# Combine both files
cat ~/bioinfo_day1/variants1.csv \
    ~/bioinfo_day1/variants2.csv \
    > ~/bioinfo_day1/all_variants.csv

echo "Combined file:"
cat ~/bioinfo_day1/all_variants.csv

# Filter for Alpha variant
echo "Alpha variant entries:"
grep "Alpha" ~/bioinfo_day1/all_variants.csv

# Count Alpha entries with a pipe
echo "Total Alpha samples:"
grep "Alpha" ~/bioinfo_day1/all_variants.csv | wc -l
```

---

## 17. Command Chaining

Run multiple commands in sequence using chaining operators.

| Operator | Behaviour |
|----------|-----------|
| `;` | Run next command regardless of success or failure |
| `&&` | Run next command **only if** previous succeeded |
| `\|\|` | Run next command **only if** previous failed |

```bash
# Semicolon: always runs both
mkdir -p ~/bioinfo_day1/test_seq ; cd ~/bioinfo_day1/test_seq ; pwd

# &&: stop if any command fails
mkdir -p ~/bioinfo_day1/success_test \
  && echo "Directory created!" \
  && echo "Ready to proceed!"

# ||: run on failure (error handling)
false || echo "Previous command failed — showing this message"
```

---

## 18. Best Practices

### Getting Help

```bash
man <command>          # Full manual pages
<command> --help       # Quick inline help
<command> -h           # Short help (tool-dependent)
```

### File & Directory Habits

```bash
# Use descriptive names
analysis_results_2024.txt    # ✅ clear
data.txt                     # ❌ vague

# Organize projects like this:
project/
├── data/
├── results/
├── scripts/
└── README.txt

# Always keep a backup before editing
cp important_file.txt important_file.bak
```

### Safety with `rm`

```bash
# Preview what you're about to delete
ls *.old

# Use interactive mode for safety
rm -i *.old
rm -r -i old_directory/
```

### Scripting & Reproducibility

```bash
# Add comments to all scripts
# Step 1: align reads to reference
bwa mem ref.fa reads.fastq > aligned.sam

# Use version control
git init
git add .
git commit -m "Add Day 1 analysis scripts"
```

---

## 19. Practice Exercises

### Exercise 1 — Directory Backup
```bash
cp -r sequencing sequencing_backup
```

### Exercise 2 — Wildcards
```bash
# Matches ethane.pdb and methane.pdb
ls *t??ne.pdb
```

### Exercise 3 — Redirection
```bash
ls run1 > sequencing_files.txt
ls run2 >> sequencing_files.txt
```

### Exercise 4 — Filter and Count
```bash
cat *_variants.csv > all_countries.csv
grep "Alpha" all_countries.csv > alpha.csv
wc -l alpha.csv
```

### Exercise 5 — Full Pipeline

```bash
echo "=== PRACTICE PROJECT ==="

# 1. Create project structure
mkdir -p ~/bioinfo_day1/practice_project/{data,results,scripts}
echo "Project structure created"

# 2. Create sample data
cat > ~/bioinfo_day1/practice_project/data/samples.csv << 'EOF'
SampleID,Variant,Count,Date
S001,Alpha,150,2024-01-01
S002,Beta,200,2024-01-02
S003,Alpha,180,2024-01-03
S004,Gamma,120,2024-01-04
S005,Alpha,190,2024-01-05
EOF
echo "Sample data created"

# 3. Filter for Alpha
grep "Alpha" ~/bioinfo_day1/practice_project/data/samples.csv \
  > ~/bioinfo_day1/practice_project/results/alpha_samples.csv
echo "Filtered for Alpha variant"

# 4. Report results
echo ""
echo "Total Alpha samples:"
wc -l ~/bioinfo_day1/practice_project/results/alpha_samples.csv

echo ""
echo "Alpha variant details:"
cat ~/bioinfo_day1/practice_project/results/alpha_samples.csv
```

---

## 20. Summary

### Commands Reference

| Category | Commands |
|----------|----------|
| Navigation | `pwd`, `ls`, `cd` |
| File management | `mkdir`, `touch`, `cp`, `mv`, `rm`, `rmdir` |
| File viewing | `cat`, `head`, `tail`, `less`, `more` |
| Searching | `grep` |
| Counting | `wc` |
| Redirection | `>`, `>>`, `<` |
| Piping | `\|` |
| Wildcards | `*`, `?` |
| Chaining | `;`, `&&`, `\|\|` |

## 📚 Commands Reference Sheet

### Navigation
```bash
pwd                    # Print working directory
ls, ls -l, ls -la      # List files (long format, show hidden)
cd <dir>, cd .., cd ~  # Change directory
```

### File Management
```bash
mkdir <dir>            # Create directory
mkdir -p <path>        # Create nested directories
touch <file>           # Create empty file
cp <src> <dest>        # Copy file
cp -r <src> <dest>     # Copy directory
mv <old> <new>         # Move/rename
rm <file>              # Delete file
rm -i <file>           # Delete with confirmation
rm -r <dir>            # Delete directory
```

### Viewing & Searching
```bash
cat <file>             # View entire file
head -n N <file>       # First N lines
tail -n N <file>       # Last N lines
less <file>            # Paginated viewer
grep <pattern> <file>  # Search for pattern
grep -v <pattern>      # Exclude pattern
grep -i <pattern>      # Case-insensitive
grep -n <pattern>      # Show line numbers
wc -l <file>           # Count lines
```

### Wildcards
```bash
*.txt                  # All .txt files
p*.pdb                 # Files starting with 'p', ending in .pdb
sample?.txt            # sample1.txt, sample2.txt, etc.
```

### Redirection & Pipes
```bash
command > file         # Write to file (overwrite)
command >> file        # Append to file
command1 | command2    # Pipe output
grep "text" file | wc -l   # Count matches
```

### Help
```bash
man <command>          # Full manual
<command> --help       # Quick help
<command> -h           # Short help
```

---

### What's Next

- Shell scripting (`if`, `for`, `while`, functions)
- Text processing (`awk`, `sed`, `cut`, `sort`, `uniq`)
- File permissions (`chmod`, `chown`)
- Process management (`ps`, `top`, `kill`)
- Remote computing (`ssh`, `scp`, `rsync`)
- Bioinformatics tools (`samtools`, `bedtools`, `bwa`, `fastqc`)

---


*Organized by: CSIR–Institute of Genomics & Integrative Biology (IGIB)*
