@echo off
echo Creating GitHub repo and pushing...
cd /d "%~dp0"

gh auth status 2>nul
if errorlevel 1 (
    echo.
    echo You need to log in to GitHub first. Run:
    echo   gh auth login -h github.com -p https -w
    echo.
    echo Then run this script again.
    pause
    exit /b 1
)

echo Creating repository drmahmoodhachim-gif.github.io...
gh repo create drmahmoodhachim-gif/drmahmoodhachim-gif.github.io --public --source=. --remote=origin --push

if errorlevel 1 (
    echo.
    echo If the repo already exists, just pushing...
    git push -u origin main
)

echo.
echo Done! Your site: https://drmahmoodhachim-gif.github.io
pause
