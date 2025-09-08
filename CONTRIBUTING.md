# Contributing to BIO559R Tutorial

Thank you for your interest in contributing to the BIO559R Tutorial! This document provides guidelines for contributing to this educational resource.

## How to Contribute

### Reporting Issues

If you encounter problems with the tutorial:

1. **Check existing issues** first to avoid duplicates
2. **Use the issue template** when creating new issues
3. **Provide detailed information**:
   - Operating system and version
   - Python/R versions
   - Error messages (full stack trace)
   - Steps to reproduce the issue

### Suggesting Improvements

We welcome suggestions for:
- Additional analysis examples
- New packages or tools to include
- Clarifications to existing documentation
- Performance improvements

### Contributing Code

1. **Fork the repository**
2. **Create a feature branch**: `git checkout -b feature/your-feature-name`
3. **Make your changes**
4. **Test your changes**:
   - Run the test script: `python scripts/test_installation.py`
   - Test any new notebooks or scripts
   - Verify documentation is accurate
5. **Commit your changes**: Use clear, descriptive commit messages
6. **Push to your fork**: `git push origin feature/your-feature-name`
7. **Submit a pull request**

### Documentation Guidelines

- Use clear, beginner-friendly language
- Include code examples where appropriate
- Test all commands and code snippets
- Follow the existing documentation style
- Include references for external resources

### Code Guidelines

- **Python**: Follow PEP 8 style guidelines
- **R**: Follow tidyverse style guide
- **Shell scripts**: Use bash and include error handling
- **Notebooks**: Include clear markdown explanations
- **Comments**: Write clear, helpful comments

### Testing

Before submitting contributions:

1. **Test installation scripts** on a clean environment
2. **Verify all notebook cells run** without errors
3. **Check documentation links** are working
4. **Run the test suite** if applicable

## Development Setup

To set up a development environment:

```bash
# Clone your fork
git clone https://github.com/your-username/BIO559R-Tutorial.git
cd BIO559R-Tutorial

# Create the conda environment
conda env create -f environment.yml
conda activate bio559r

# Test the installation
python scripts/test_installation.py
```

## Pull Request Process

1. **Update documentation** if you change functionality
2. **Add tests** for new features when applicable
3. **Update the changelog** if significant changes are made
4. **Ensure CI passes** (when implemented)
5. **Request review** from maintainers

## Code of Conduct

This project follows a simple code of conduct:

- **Be respectful** and inclusive
- **Be constructive** in feedback
- **Focus on education** and helping students learn
- **Acknowledge contributions** from others

## Questions?

If you have questions about contributing:

- Open an issue with the "question" label
- Check existing documentation first
- Be specific about what you need help with

## Recognition

Contributors will be acknowledged in:
- The main README file
- Release notes for significant contributions
- The contributors section (when implemented)

Thank you for helping make this tutorial better for all students! 🎓

