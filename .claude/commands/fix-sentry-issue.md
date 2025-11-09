Please analyze and fix Sentry issueId: $ARGUMENTS.

organizationSlug: talley-lambert
projectSlugOrId: fpbase

Follow these steps:

1. Use sentry mcp with `get_issue_details` and maybe `get_trace_details` to get information about the issue
2. Understand the problem described in the issue
3. Search the codebase for relevant files
4. Plan the necessary changes to fix the issue

Stop here, and present your understanding of the issue and your proposed fix before proceeding.

5. Implement the necessary changes to fix the issue
6. Write and run tests to verify the fix
7. Ensure code passes linting and type checking
8. Stop there, don't commit or push changes.
9. Use sentry mcp with `search_issues` to find any similar issues that this fix might also close

CRITICAL!
Do NOT just silence the error; you MUST understand and fix the root cause!
If you don't understand the issue, use web search to research it.
If you still don't know why the issue is happening, stop there.
