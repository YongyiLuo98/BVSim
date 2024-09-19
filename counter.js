let visitCount = localStorage.getItem('visitCount') || 0;
visitCount++;
localStorage.setItem('visitCount', visitCount);
document.getElementById('visitCount').innerText = `访问次数: ${visitCount}`;
