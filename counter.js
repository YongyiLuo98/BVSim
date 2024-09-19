let visitCount = localStorage.getItem('visitCount') || 0;
visitCount++;
localStorage.setItem('visitCount', visitCount);
document.getElementById('visitCount').innerText = `BVSim profile visits: ${visitCount}`;
