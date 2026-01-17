# VS Code Debug Configuration for VCF Tests

This directory contains VS Code launch configurations for debugging the VCF transformer tests.

## Available Debug Configurations

### 1. Debug All Tests
Runs all 12 transformer tests with debugger attached.
- **Usage**: Select "Debug All Tests" from the debug dropdown
- **Command equivalent**: `python -m unittest publisher.VCF_test_transformers -v`

### 2. Debug Focused Test (FOCUS_TEST=true)
Runs only the test marked with `@focus` decorator (currently TranscriptsTransformer).
- **Usage**: Select "Debug Focused Test (FOCUS_TEST=true)" from the debug dropdown
- **Command equivalent**: `FOCUS_TEST=true python -m publisher.VCF_test_transformers`
- **Note**: Move the `@focus` decorator in `VCF_test_transformers.py` to focus on different tests

### 3. Debug Specific Transformer Tests
Individual configurations for debugging specific transformer tests:
- Debug Genes Transformer Test
- Debug Transcripts Transformer Test
- Debug Variants Transformer Test

Add more as needed following the pattern in `launch.json`.

### 4. Debug Current Test File
Generic configuration to debug whatever Python file is currently open.
- **Usage**: Open the test file and select "Debug Current Test File"

## How to Use

1. **Set breakpoints**: Click in the gutter left of line numbers in your Python files
2. **Select configuration**: Use the dropdown in the Debug view (Ctrl+Shift+D / Cmd+Shift+D)
3. **Start debugging**: Press F5 or click the green play button
4. **Debug controls**:
   - Continue (F5)
   - Step Over (F10)
   - Step Into (F11)
   - Step Out (Shift+F11)
   - Restart (Ctrl+Shift+F5)
   - Stop (Shift+F5)

## Important Notes

- **No import changes needed**: These configurations use the `-m` module syntax, so imports like `import publisher.VCF_transformers` work correctly
- **Working directory**: Set to `${workspaceFolder}` (the repository root)
- **justMyCode**: Set to `false` so you can step into library code if needed
- **Python environment**: Make sure you have the correct Python interpreter selected in VS Code

## Testing Without Debugger

To verify the configurations work, you can run the equivalent commands in terminal:

```bash
# Run all tests
python -m unittest publisher.VCF_test_transformers -v

# Run focused test
FOCUS_TEST=true python -m publisher.VCF_test_transformers

# Run specific test class
python -m unittest publisher.VCF_test_transformers.TestTranscriptsTransformer -v
```

## Troubleshooting

If imports fail:
1. Make sure your Python interpreter is set correctly in VS Code
2. Check that `python.analysis.extraPaths` includes `${workspaceFolder}` in `.vscode/settings.json`
3. Verify the working directory is set to `${workspaceFolder}` in the launch configuration
