# Claude Code Setup for FPbase

This directory contains configuration for running FPbase in Claude Code's cloud environment.

## Automatic Setup

The setup script runs automatically on session start via hooks. It configures:

- ✅ PostgreSQL database
- ✅ Python dependencies (uv)
- ✅ Node dependencies (pnpm)
- ✅ Heroku CLI
- ✅ MCP Servers (GitHub, Sentry, Playwright)

## Required Environment Variables

To enable all features, add these environment variables to your Claude Code project settings:

### GitHub MCP Server
```bash
GITHUB_TOKEN=ghp_your_personal_access_token_here
```
**Where to get it:**
- Go to https://github.com/settings/tokens
- Generate a personal access token with `repo` scope
- **Capabilities:** Create/manage issues, PRs, view repo info, manage releases

### Sentry MCP Server
```bash
SENTRY_AUTH_TOKEN=your_sentry_auth_token_here
SENTRY_ORG=your_sentry_org_slug
```
**Where to get it:**
- Sentry Auth Token: https://sentry.io/settings/account/api/auth-tokens/
- Create a token with `project:read`, `event:read`, `org:read` scopes
- Sentry Org: Your organization slug (e.g., `fpbase`)
- **Capabilities:** View errors, search issues, get project stats

### Heroku CLI
```bash
HEROKU_API_KEY=your_heroku_api_key_here
```
**Where to get it:**
- Go to https://dashboard.heroku.com/account
- Scroll to "API Key" section
- Click "Reveal" to see your API key
- **Capabilities:** Deploy apps, manage dynos, view logs, run commands

## MCP Servers Configured

### 1. GitHub (`@modelcontextprotocol/server-github`)
- Create and manage issues
- Create and manage pull requests
- Search code and repositories
- List and manage branches
- View repository information

### 2. Sentry (`@modelcontextprotocol/server-sentry`)
- Search and view error issues
- Get issue details and statistics
- List projects and their health
- View error trends

### 3. Playwright (`@executeautomation/playwright-mcp-server`)
- Run browser automation scripts
- Execute e2e tests locally
- Take screenshots and videos
- Test responsive designs

## Manual Setup (if needed)

If the automatic setup fails, run manually:

```bash
bash .claude/scripts/setup_cloud_environment.sh
```

## Testing

Run tests with:

```bash
just test-py   # Python tests (165 tests)
just test-js   # JavaScript tests (3 tests)
just test-e2e  # End-to-end tests (49 tests)
just test      # All tests
```

## Heroku CLI Usage

Once `HEROKU_API_KEY` is set, you can use Heroku CLI:

```bash
heroku login -i  # Login with API key
heroku apps      # List your apps
heroku logs -a fpbase --tail  # View logs
heroku run bash -a fpbase     # Run commands on production
```

## Troubleshooting

### PostgreSQL not starting
The setup script fixes common PostgreSQL issues automatically:
- SSL configuration
- Authentication methods
- File permissions

### MCP servers not working
Check that environment variables are set:
```bash
echo $GITHUB_TOKEN
echo $SENTRY_AUTH_TOKEN
echo $HEROKU_API_KEY
```

### Heroku CLI authentication
If you get auth errors:
```bash
export HEROKU_API_KEY=your_key_here
heroku auth:whoami
```

## Architecture

```
.claude/
├── settings.json              # MCP servers and hooks config
├── scripts/
│   └── setup_cloud_environment.sh  # Automatic environment setup
└── README-CLAUDE-SETUP.md     # This file
```

## Security Note

Never commit API keys or tokens to the repository. They should only be set as environment variables in your Claude Code project settings.
